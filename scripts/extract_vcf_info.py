"""
extract_vcf_info.py

This script extracts specific information from a VCF (Variant Call Format) file and writes the output to a CSV file. 
The information that can be extracted includes SNP information (SNP_INFO), Genotype (GT), Estimated Alternate Allele Dosage (DS), or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP).

Usage:
    python3 extract_vcf_info.py --vcf <vcf_file> --extract <info_type> --out <output_file>

Arguments:
    --vcf: The VCF file to parse. This is a required argument.
    --extract: The type of information to extract. Choices are comma separated: SNP_INFO, GT, DS, GP. If not specified, it will only generate the SNP_INFO file.
    --out: The name of the output file (including the file path). This is a required argument.

Example:
    python3 extract_vcf_info.py --vcf chr1.vcf.gz --extract GT,SNP_INFO --out chr1_info

This will extract the Genotype (GT) and SNP information (SNP_INFO) from the chr1.vcf.gz file and write the output to chr1_info_GT.csv and chr1_info_SNP_INFO.csv respectively.

Note: Ensure that the VCF file exists and the output file can be written.
"""

import pysam
import argparse
import os
import pandas as pd

def extract_vcf_info(vcf_file, extract):
    """
    Extracts SNP information such as chromosome, position, reference allele, alternate allele, allele frequency (AF), minor allele frequency (MAF), imputation accuracy (R2), empirical R-square (ER2), and additional information including (IMPUTED, TYPED, TYPED_ONLY), Genotype (GT), Estimated Alternate Allele Dosage (DS) or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP) information from a VCF file and writes it to a CSV file.

    Args:
        vcf_file (str): Path to the VCF file.
        extract (list): A list containing any of the following: ["SNP_INFO", "GT","DS","GP"]

    Returns:
        None
    """
    print(f"Parsing file: {vcf_file}")
    for option in extract:
        output_file_opt = f"{args.out}_{option}.tsv"
        with pysam.VariantFile(vcf_file) as vcf:
            with open(output_file_opt, 'w') as out:
                sample_names = vcf.header.samples
                headers = ['SNP','chr', 'position', 'REF', 'ALT','AF','MAF','Rsq','RSqEmp','INFO']
                if option!= 'SNP_INFO':
                    out.write("IID\t" + '\t'.join(sample_names) + "\n")
                else:
                    out.write('\t'.join(headers) + "\n")
                for record in vcf:
                    chrom = "NA"
                    pos = "NA"
                    ref = "NA"
                    alt = "NA"
                    info = []
                    extra_info = "NA"
                    if record.chrom and record.pos and record.ref and record.alts:
                        chrom = record.chrom
                        pos = record.pos
                        ref = record.ref
                        alt = ','.join(map(str, record.alts))
                    if option == 'SNP_INFO':
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
                    if option == 'GT':
                        info = ['|'.join(map(str, sample['GT'])) if sample['GT'] is not None else 'NA' for sample in record.samples.values()]
                        out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")
                    elif option == 'DS':
                        info = [str(int(sample['DS'])) if sample['DS'].is_integer() else '{:.3f}'.format(sample['DS']) if sample['DS'] is not None else 'NA' for sample in record.samples.values()]
                        out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")
                    elif option == 'GP':
                        info = ['|'.join(str(int(value)) if value.is_integer() else f'{value:.3f}' for value in sample['GP']) if sample['GP'] is not None else 'NA' for sample in record.samples.values()]
                        out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")
        if option!= 'SNP_INFO':
            # transposed file:
            df = pd.read_csv(output_file_opt, sep='\t', index_col=0)
            df_transposed = df.transpose().reset_index()
            df_transposed.rename(columns={'index': 'IID'}, inplace=True)
            df_transposed.to_csv(f"{args.out}_{option}.csv", sep=',', index=False)
            os.remove(output_file_opt)
            print(f"Output generated: {args.out}_{option}.csv")
        else:
            print(f"Output generated: {output_file_opt}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='This script can extract SNP information (SNP_INFO), Genotype (GT), Estimated Alternate Allele Dosage (DS), or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP) information from a VCF file and writes the output to a CSV file.',
        epilog='Ensure that the VCF file exists and the output file can be written.\n\nExample usage:\npython3 extract_vcf_info.py --vcf chr1.vcf.gz --extract GT,SNP_INFO --out chr1_info\n',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--vcf', help='The VCF file to parse.', required=True)
    parser.add_argument('--extract', help='The type of information to extract. Choices are comma separated: SNP_INFO, GT, DS, GP. If not specified, it will only generate the SNP_INFO file.', default=None)
    parser.add_argument('--out', help='Provide the base name of output files (including the file path). The script will append the type of information to the base name.', required=True)
    args = parser.parse_args()

    valid_options = ['SNP_INFO', 'GT', 'DS', 'GP']
    if args.extract:
        user_options = args.extract.split(',')
        extract = ["SNP_INFO"] + [option for option in user_options if option in valid_options and option != 'SNP_INFO']
    else:
        extract = ["SNP_INFO"]

    extract_vcf_info(args.vcf, extract)