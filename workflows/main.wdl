version 1.0

import "./tasks/extract_info.wdl" as info

workflow main {
    
    input {
        File query_variants
        File? query_samples
        Array [File] imputed_vcf
        String prefix
    }
    
    parameter_meta {
        query_variants: "A tab delimited file containing the list of query variants. Each line should have the format: Chromosome, Pos, ID, Ref, Alt."
        query_samples: "A file containing a list of sample IDs. The file should have one sample ID per row. (optional)"
        imputed_vcf: "An array of imputed VCF files and their indices. The VCF files should be in .vcf.gz format and the indices should be in CSI or TBI format."
        prefix: "A prefix for the output files."
    }

    call info.extract {
        input: query_variants = query_variants, query_samples = query_samples, imputed_vcf = imputed_vcf, prefix = prefix
    }
    
    output {
        File snp_info = extract.snp_info
        File genotype_info = extract.gt_info
        File dosage_info  = extract.ds_info
        File geno_prob_info = extract.gp_info
    }
    
    meta {
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
        description: "This workflow extracts information from a set of imputed VCF files using a list of query variants or sample IDs. It can extract the following information: chromosome, position, reference allele, alternate allele, allele frequency (AF), minor allele frequency (MAF), imputation accuracy (R2), empirical R-square (ER2), Genotype (GT), Estimated Alternate Allele Dosage (DS), and Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP). The output is a set of files containing the extracted information."
    }
}