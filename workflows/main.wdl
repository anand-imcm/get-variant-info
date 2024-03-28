version 1.0

import "./tasks/extract_info.wdl" as info

workflow main {
    
    String pipeline_version = "1.0.1"
    String container_src = "ghcr.io/anand-imcm/get-variant-info:~{pipeline_version}"
    
    input {
        File query_variants
        File? query_samples
        Array [File] imputed_vcf
        String prefix
        String extract_item = "GT"
    }
    
    parameter_meta {
        query_variants: "A tab delimited file with list of query variants. Each line should be formatted as: Chromosome, Pos, ID, Ref, Alt. (required)"
        query_samples: "A file with list of sample IDs. Each line should contain one sample ID. (optional)"
        imputed_vcf: "Array of imputed VCF files and their indices. VCF files should be in .vcf.gz format and indices in CSI or TBI format. (required)"
        prefix: "Prefix for the output files. (required)"
        extract_item: "A string specifying the information to extract from the FORMAT field of the VCF file. The available choices are GT, DS, and GP. Please provide as a comma-separated string. If not provided, it defaults to 'GT'."
    }
    
    call info.extract {
        input: query_variants = query_variants, query_samples = query_samples, imputed_vcf = imputed_vcf, prefix = prefix, extract_item = extract_item, docker = container_src
    }
    
    output {
        File snp_info = extract.snp_info
        File? genotype_info = extract.gt_info
        File? dosage_info  = extract.ds_info
        File? geno_prob_info = extract.gp_info
    }
    
    meta {
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
        description: "This workflow extracts information from a set of imputed VCF files using a list of query variants or sample IDs. It can extract the following information: chromosome, position, reference allele, alternate allele, allele frequency (AF), minor allele frequency (MAF), imputation accuracy (R2), empirical R-square (ER2), Genotype (GT), Estimated Alternate Allele Dosage (DS), and Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP). The output is a set of files containing the extracted information."
    }
}