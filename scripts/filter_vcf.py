"""
filter_vcf.py

This script filters a VCF file based on a TSV file and writes the filtered variants to a new VCF file.

Usage:
    python filter_vcf.py --query <path_to_tsv_file> --vcf <path_to_vcf_file> --out <base_name_of_output_file>

Example:
    python filter_vcf.py --query data.tsv --vcf variants.vcf --out filtered

This will create a new VCF file named "filtered.vcf" that only contains the variants that are present in both "variants.vcf" and "data.tsv".

The TSV file should be tab-separated and have the following columns (without headers):
    1. Chromosome
    2. Position
    3. ID
    4. Reference allele
    5. Alternate alleles (comma-separated if there are multiple)

The VCF file should be in standard VCF format.
"""

import pandas as pd
import pysam
import argparse

def filter_vcf(tsv_file, vcf_file, out_prefix):
    """
    Filters a VCF file based on a TSV file and writes the filtered variants to a new VCF file.

    Parameters:
    tsv_file (str): The path to the TSV file to use for filtering.
    vcf_file (str): The path to the VCF file to filter.
    out_prefix (str): The base name of the output file.
    """
    df = pd.read_csv(tsv_file, sep='\t', names=['chr', 'pos', 'id', 'ref', 'alt'])
    vcf_in = pysam.VariantFile(vcf_file)
    vcf_out = pysam.VariantFile(f"{out_prefix}.vcf", 'w', header=vcf_in.header)
    for record in vcf_in:
        matching_rows = df[(df['chr'] == record.chrom) & (df['pos'] == record.pos) & (df['ref'] == record.ref)]
        for _, row in matching_rows.iterrows():
            for alt in record.alts:
                if str(alt) in row['alt'].split(','):
                    new_record = record.copy()
                    new_record.alts = [alt]
                    vcf_out.write(new_record)
                    break
    vcf_out.close()
    vcf_in.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter a VCF file based on a TSV file.')
    parser.add_argument('--query', type=str, help='(Required) A tab delimited file with list of query variants. Each line should be formatted as: Chromosome, Pos, ID, Ref, Alt.')
    parser.add_argument('--vcf', type=str, help='(Required) The VCF file to filter.')
    parser.add_argument('--out', type=str, help='(Required) The base name of the output file.')
    args = parser.parse_args()

    filter_vcf(args.query, args.vcf, args.out)