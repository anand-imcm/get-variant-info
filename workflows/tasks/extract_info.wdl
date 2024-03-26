version 1.0

task extract_info {
    input {
        File query_variants
        File? query_samples
        File imputed_vcf
        String prefix
    }

    command <<<
        if [ -f ~{query_variants} ] && [ -f ~{query_samples} ]; then
            bcftools view --regions-file ~{query_variants} --samples-file ~{query_samples} ~{imputed_vcf} > ~{prefix}_query_subset.vcf
        else
            bcftools view --regions-file ~{query_variants} ~{imputed_vcf} > ~{prefix}_query_subset.vcf
        fi
        # extract snp INFO
        python3 extract_snp_info.py --vcf ~{prefix}_query_subset.vcf --out ~{prefix}_query_subset_extracted_snps_info.tsv
        # extract FORMAT fields
        python3 /home/anand/Documents/aspire-files/data-oxford/terra.bio/get-variant-info/scripts/extract_vcf_info.py --vcf ~{prefix}_query_subset.vcf --out ~{prefix}_extracted_snps_bg_genotype.csv --extract GT
        python3 /home/anand/Documents/aspire-files/data-oxford/terra.bio/get-variant-info/scripts/extract_vcf_info.py --vcf ~{prefix}_query_subset.vcf --out ~{prefix}_extracted_snps_dosage.csv --extract DS
        python3 /home/anand/Documents/aspire-files/data-oxford/terra.bio/get-variant-info/scripts/extract_vcf_info.py --vcf ~{prefix}_query_subset.vcf --out ~{prefix}_extracted_snps_bg_prob.csv --extract GP
    >>>
    output {
        File snp_info = prefix + "_query_subset_extracted_snps_info.tsv"
        File gt_info = prefix + "_extracted_snps_bg_genotype.csv"
        File ds_info = prefix + "_extracted_snps_dosage.csv"
        File gp_info = prefix + "_extracted_snps_bg_prob.csv"
    }
}