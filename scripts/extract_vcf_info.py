"""
extract_vcf_info.py

This script extracts specific information from a VCF (Variant Call Format) file and writes the output to a CSV file. 
The information that can be extracted includes SNP information (SNP_INFO), Genotype (GT), Estimated Alternate Allele Dosage (DS), or Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP).

Usage:
    python3 extract_vcf_info.py --vcf <vcf_file> --extract <info_type> --out <output_file> [--ped <ped_file>]

Arguments:
    --vars: The list of query variants. Each line should be formatted as: Chromosome, Pos, ID, Ref, Alt. This will be used to filter the output.
    --vcf: The VCF file to parse. This is a required argument.
    --extract: The type of information to extract. Choices are comma separated: SNP_INFO, GT, DS, GP. If not specified, it will only generate the SNP_INFO file.
    --out: The name of the output file (including the file path). This is a required argument.
    --ped: (Optional) A .ped file containing compound genotypes. If supplied, genotypes will be sourced from this file. Otherwise, genotypes will be extracted from the VCF file.

Example:
    python3 extract_vcf_info.py --vcf chr1.vcf.gz --extract GT,SNP_INFO --out chr1_info --ped chr1.ped

This will extract the Genotype (GT) and SNP information (SNP_INFO) from the chr1.vcf.gz file and write the output to chr1_info_GT.csv and chr1_info_SNP_INFO.csv respectively. If a .ped file is provided, the genotypes will be sourced from this file.

Note: Ensure that the VCF file exists and the output file can be written.
"""

import pysam
import argparse
import os
import pandas as pd

def get_query_vars(query_list):
    """
    Extracts chr:pos:ref:alt from the given file and returns a unique list of variants.

    Args:
        query_list (str) : The list of query variants. Each line should be formatted as: Chromosome, Pos, ID, Ref, Alt. This will be used to filter the output

    Returns:
        unique_ids (list) :  A list containing variants in chr:pos:ref:alt format.
    """
    query_df = pd.read_csv(query_list, sep='\t', header=None, names=['chr', 'pos', 'source', 'ref', 'alt'])
    query_df['id'] = query_df['chr'].astype(str) + ":" + query_df['pos'].astype(str) + ":" + query_df['ref'] + ":" + query_df['alt']
    unique_ids = query_df['id'].unique().tolist()
    return unique_ids

def extract_vcf_info(vcf_file, extract, ped, query_variants):
    """
    Extracts specific information from a VCF file and writes it to a CSV file. 
    The information that can be extracted includes SNP information (SNP_INFO), 
    Genotype (GT), Estimated Alternate Allele Dosage (DS), or Estimated Posterior 
    Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP). If a .ped file is provided, 
    genotypes will be sourced from this file. Otherwise, genotypes will be extracted 
    from the VCF file.

    Args:
        vcf_file (str): Path to the VCF file.
        extract (list): A list containing any of the following: ["SNP_INFO", "GT","DS","GP"]
        ped (str): Optional; Path to the .ped file containing compound genotypes.
        query_variants (list) : A list containing variants in chr:pos:ref:alt format.

    Returns:
        None
    """
    print(f"Parsing file: {vcf_file}")
    for option in extract:
        output_file_opt = f"{args.out}_{option}.tsv"
        all_variant_IDs = []
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
                        vcf_var = f"{chrom}:{pos}:{ref}:{alt}"
                        if vcf_var in query_variants:
                            all_variant_IDs.append(f"{chrom}:{pos}:{ref}:{alt}")
                        
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
                        if f"{chrom}:{pos}:{ref}:{alt}" in query_variants:
                            out.write(f"{chrom}:{pos}:{ref}:{alt}\t{chrom}\t{pos}\t{ref}\t{alt}\t{af}\t{maf}\t{r2}\t{er2}\t{extra_info}\n")
                    if option == 'GT' and not ped:
                        info = []
                        for sample in record.samples.values():
                            if sample['GT'] is not None:
                                if sample['GT'] == (0, 0):
                                    info.append(ref + ref)
                                elif sample['GT'] in [(0, 1), (1, 0)]:
                                    info.append(ref + alt)
                                elif sample['GT'] == (1, 1):
                                    info.append(alt + alt)
                            else:
                                info.append('NA')
                        if f"{chrom}:{pos}:{ref}:{alt}" in query_variants:
                            out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")
                    if option == 'DS':
                        info = [str(int(sample['DS'])) if sample['DS'].is_integer() else '{:.3f}'.format(sample['DS']) if sample['DS'] is not None else 'NA' for sample in record.samples.values()]
                        if f"{chrom}:{pos}:{ref}:{alt}" in query_variants:
                            out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")
                    if option == 'GP':
                        info = ['|'.join(str(int(value)) if value.is_integer() else f'{value:.3f}' for value in sample['GP']) if sample['GP'] is not None else 'NA' for sample in record.samples.values()]
                        if f"{chrom}:{pos}:{ref}:{alt}" in query_variants:
                            out.write(f"{chrom}:{pos}:{ref}:{alt}\t" + '\t'.join(info) + "\n")
            ped_headers = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'Phenotype'] + all_variant_IDs
            ped_columns_to_ignore = ['FID', 'PAT', 'MAT', 'SEX', 'Phenotype']
            if option == 'GT' and ped:
                ped_df_all = pd.read_csv(ped, sep='\t', header=None, names=ped_headers)
                ped_df_all = ped_df_all.drop(columns=ped_columns_to_ignore)
                ped_df_all.to_csv(f"{args.out}_{option}.csv", sep=',', index=False)
                print(f"Output generated: {args.out}_{option}.csv")
                os.remove(output_file_opt)
            else:
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
    parser.add_argument('--vars', help='(Optional) List of query variants. Each line should be formatted as: Chromosome, Pos, ID, Ref, Alt. This will be used to filter the output.',default=None)
    parser.add_argument('--vcf', help='(Required) The VCF file to parse.', required=True)
    parser.add_argument('--extract', help='(Required) The type of information to extract. Choices are comma separated: SNP_INFO, GT, DS, GP. If not specified, it will only generate the SNP_INFO file.', default=None)
    parser.add_argument('--ped', help='(Optional) Provide a .ped file containing compound genotypes. If supplied, genotypes will be sourced from this file. Otherwise, genotypes will be extracted from the VCF file.', default=None)
    parser.add_argument('--out', help='(Required) Provide the base name of output files (including the file path). The script will append the type of information to the base name.', required=True)
    args = parser.parse_args()

    valid_options = ['SNP_INFO', 'GT', 'DS', 'GP']
    if args.extract:
        user_options = args.extract.split(',')
        extract = ["SNP_INFO"] + [option for option in user_options if option in valid_options and option != 'SNP_INFO']
    else:
        extract = ["SNP_INFO"]

    query_vars = get_query_vars(args.vars)
    extract_vcf_info(args.vcf, extract, args.ped, query_vars)