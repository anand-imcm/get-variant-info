"""
extract_vcf_info.py

This script extracts specific information from a VCF (Variant Call Format) file and writes the output to a CSV file. 
The information that can be extracted includes Genotype (GT), Estimated Alternate Allele Dosage (DS), or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP).

Usage:
    python3 extract_vcf_info.py --vcf <vcf_file> --extract <info_type> --out <output_file>

Arguments:
    --vcf: The VCF file to parse. This is a required argument.
    --extract: The type of information to extract. Choices are 'GT', 'DS', 'GP'. This is a required argument.
    --out: The name of the output file (including the file path). This is a required argument.

Example:
    python3 extract_vcf_info.py --vcf chr1.vcf.gz --extract GT --out chr1_info.csv

This will extract the Genotype (GT) information from the chr1.vcf.gz file and write the output to chr1_info.csv.

Note: Ensure that the VCF file exists and the output file can be written.
"""

import pysam
import argparse
import os
import pandas as pd

def extract_vcf_info(vcf_file, extract, output_file):
    """
    Extracts Genotype (GT), Estimated Alternate Allele Dosage (DS) or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP) information from a VCF file and writes it to a CSV file.

    Args:
        vcf_file (str): Path to the VCF file.
        output_file (str): Path to the output file.

    Returns:
        None
    """
    print(f"Parsing file: {vcf_file}")
    with pysam.VariantFile(vcf_file) as vcf:
        with open(output_file, 'w') as out:
            sample_names = vcf.header.samples
            out.write("IID\t" + '\t'.join(sample_names) + "\n")
            for record in vcf:
                chrom = "NA"
                pos = "NA"
                ref = "NA"
                alt = "NA"
                info = []
                if record.chrom and record.pos and record.ref and record.alts:
                    chrom = record.chrom
                    pos = record.pos
                    ref = record.ref
                    alt = ','.join(map(str, record.alts))
                if extract == 'GT':
                    info = ['|'.join(map(str, sample['GT'])) if sample['GT'] is not None else 'NA' for sample in record.samples.values()]
                elif extract == 'DS':
                    info = [str(int(sample['DS'])) if sample['DS'].is_integer() else '{:.3f}'.format(sample['DS']) if sample['DS'] is not None else 'NA' for sample in record.samples.values()]
                elif extract == 'GP':
                    info = ['|'.join(str(int(value)) if value.is_integer() else f'{value:.3f}' for value in sample['GP']) if sample['GP'] is not None else 'NA' for sample in record.samples.values()]
                out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script can extract Genotype (GT) or Estimated Alternate Allele Dosage (DS) or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP) information from a VCF file writes the output to a CSV file.',
        epilog='Ensure that the VCF file exists and the output file can be written.\n\nExample usage:\npython3 extract_vcf_info.py --vcf chr1.vcf.gz --extract GT --out chr1_info.csv\n',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--vcf', help='The VCF file to parse.', required=True)
    parser.add_argument('--extract', help='The type of information to extract. Choices are "GT", "DS", "GP".', required=True, choices=['GT', 'DS', 'GP'], )
    parser.add_argument('--out', help='Provide the name of output file (including the file path)', required=True)
    args = parser.parse_args()
    output_file_temp = f"{args.out}.temp.txt"
    extract_vcf_info(args.vcf, args.extract, output_file_temp)

    # transposed file:
    df = pd.read_csv(output_file_temp, sep='\t', index_col=0)
    df_transposed = df.transpose().reset_index()
    df_transposed.rename(columns={'index': 'IID'}, inplace=True)
    df_transposed.to_csv(args.out, sep=',', index=False)
    os.remove(output_file_temp)
    print(f"Output generated: {args.out}")