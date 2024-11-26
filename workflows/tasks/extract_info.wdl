version 1.0

task extract {

    input {
        File query_variants
        File? query_samples
        Array [File] imputed_vcf
        String prefix
        String extract_item
        String docker
        Boolean use_GT_from_PED = false
        Boolean match_pos_only = false
    }
    
    Int disk_size_gb = ceil(size(imputed_vcf, "GiB")) + 5
    
    command <<<
        for vcf in ~{sep=' ' imputed_vcf}; do
            ln -s $vcf $(basename $vcf)
            if [[ $vcf == *.vcf.gz ]]; then
                if [ ! -f $vcf.csi ] && [ ! -f $vcf.tbi ]; then
                    bcftools index -c $(basename $vcf)
                fi
            fi
        done
        
        CHROMOSOMES=$(awk '{print $1}' ~{query_variants} | sort | uniq)
        
        VCF_FILES=""
        for chr in $CHROMOSOMES; do
            awk -v chr=$chr '$1 == chr' ~{query_variants} > ${chr}_query_subset_variants_list.txt
            if [ -f ${chr}_query_subset_variants_list.txt ]; then
                if [ -n "~{query_samples}" ] && [ -f ~{query_samples} ]; then
                    bcftools view --regions-file ${chr}_query_subset_variants_list.txt --samples-file ~{query_samples} --force-samples ${chr}.dose*vcf.gz > ~{prefix}_${chr}_query_subset.vcf
                else
                    bcftools view --regions-file ${chr}_query_subset_variants_list.txt ${chr}.dose*vcf.gz > ~{prefix}_${chr}_query_subset.vcf
                fi
                if [ "~{match_pos_only}" = "false" ]; then
                    python3 /scripts/filter_vcf.py --query ${chr}_query_subset_variants_list.txt --vcf ~{prefix}_${chr}_query_subset.vcf --out ~{prefix}_${chr}_query_matched_alleles_subset
                    VCF_FILES="${VCF_FILES} ~{prefix}_${chr}_query_matched_alleles_subset.vcf"
                else
                    VCF_FILES="${VCF_FILES} ~{prefix}_${chr}_query_subset.vcf"
                fi
            fi
        done
        
        bcftools concat ${VCF_FILES} -o ~{prefix}_query_extracted_concat.vcf
        bcftools norm --rm-dup exact ~{prefix}_query_extracted_concat.vcf -o ~{prefix}_query_extracted.vcf
        
        if [ "~{use_GT_from_PED}" = "true" ]; then
            plink2 --vcf ~{prefix}_query_extracted.vcf --recode compound-genotypes --out ~{prefix}_query_extracted
            python3 /scripts/extract_vcf_info.py --vcf ~{prefix}_query_extracted.vcf --out ~{prefix}_extracted --extract ~{extract_item} --ped ~{prefix}_query_extracted.ped
        else
            python3 /scripts/extract_vcf_info.py --vcf ~{prefix}_query_extracted.vcf --out ~{prefix}_extracted --extract ~{extract_item}
        fi
    >>>
    
    output {
        File vcf = prefix + "_query_extracted.vcf"
        File snp_info = prefix + "_extracted_SNP_INFO.tsv"
        File? gt_info = prefix + "_extracted_GT.csv"
        File? ds_info = prefix + "_extracted_DS.csv"
        File? gp_info = prefix + "_extracted_GP.csv"
    }
    
    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}