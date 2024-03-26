"""
extract_snp_info.py

This script parses a VCF file and extracts the chromosome, position, reference allele, alternate allele, allele frequency (AF), minor allele frequency (MAF), imputation accuracy (R2), empirical R-square (ER2), and additional information (IMPUTED, TYPED, TYPED_ONLY). The output is written to a TSV file.

Usage:
    python3 extract_snp_info.py --vcf <vcf_file> --out <output_file>

Arguments:
    --vcf: The VCF file to parse. This is a required argument.
    --out: The name of the output file (including the file path). This is a required argument.

Example:
    python3 extract_snp_info.py --vcf chr1.vcf.gz --out chr1_snp_info.tsv

This will parse the chr1.vcf.gz file, extract the SNP information, and write the output to chr1_snp_info.tsv.

Note: Ensure that the VCF file exists and the output file can be written.
"""

import pysam
import argparse

def extract_snp_info(vcf_file, output_file):
    """
    Extracts SNP information such as chromosome, position, reference allele, alternate allele, allele frequency (AF), minor allele frequency (MAF), imputation accuracy (R2), empirical R-square (ER2), and additional information (IMPUTED, TYPED, TYPED_ONLY) from a VCF file and writes it to a TSV file.

    Args:
        vcf_file (str): Path to the VCF file.
        output_file (str): Path to the output file.

    Returns:
        None
    """
    print(f"Parsing file: {vcf_file}")
    with pysam.VariantFile(vcf_file) as vcf:
        with open(output_file, 'w') as out:
            headers = ['SNP','chr', 'position', 'REF', 'ALT','AF','MAF','Rsq','RSqEmp','INFO']
            out.write('\t'.join(headers) + "\n")
            for record in vcf:
                chrom = "NA"
                pos = "NA"
                ref = "NA"
                alt = "NA"
                extra_info = "NA"
                if record.chrom and record.pos and record.ref and record.alts:
                    chrom = record.chrom
                    pos = record.pos
                    ref = record.ref
                    alt = ','.join(map(str, record.alts))
                af = 'NA' if record.info.get('AF', 'NA') == 'NA' else format(record.info.get('AF', 'NA'), '.5f')
                maf = 'NA' if record.info.get('MAF', 'NA') == 'NA' else format(record.info.get('MAF', 'NA'), '.5f')
                r2 = 'NA' if record.info.get('R2', 'NA') == 'NA' else format(record.info.get('R2', 'NA'), '.5f')
                er2 = 'NA' if record.info.get('ER2', 'NA') == 'NA' else format(record.info.get('ER2', 'NA'), '.5f')
                if 'IMPUTED' in record.info:
                    extra_info = "IMPUTED"
                if 'TYPED' in record.info:
                    extra_info = "TYPED"
                if 'TYPED_ONLY' in record.info:
                    extra_info = "TYPED_ONLY"
                out.write(f"{chrom}:{pos}:{ref}:{alt}\t{chrom}\t{pos}\t{ref}\t{alt}\t{af}\t{maf}\t{r2}\t{er2}\t{extra_info}\n")
            print(f"Output generated: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script extracts the chromosome, position, reference allele, alternate allele, allele frequency (AF), minor allele frequency (MAF), imputation accuracy (R2), empirical R-square (ER2), and additional information (IMPUTED, TYPED, TYPED_ONLY) from a VCF file and writes the output to a TSV file.',
        epilog='Ensure that the VCF file exists and the output file can be written.\n\nExample usage: python3 extract_snp_info.py --vcf chr1.vcf.gz --out chr1_snp_info.tsv',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--vcf', help='The VCF file to parse.', required=True)
    parser.add_argument('--out', help='Provide the name of output file (including the file path)', required=True)
    args = parser.parse_args()
    output_file = f"{args.out}"
    extract_snp_info(args.vcf, output_file)