version 1.0

import "./tasks/extract_info.wdl" as info

workflow main {
    input {
        File query_variants
        File? query_samples
        Array [File] imputed_vcf
        String prefix
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
}