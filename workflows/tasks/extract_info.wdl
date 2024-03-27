version 1.0

task extract {

    input {
        File query_variants
        File? query_samples
        Array [File] imputed_vcf
        String prefix
        String extract_item
    }
    
    Int disk_size_gb = ceil(size(imputed_vcf, "GiB")) + 5
    
    command <<<
        for vcf in ~{sep=' ' imputed_vcf}; do
            ln -s $vcf $(basename $vcf)
            # If the index file does not exist, create it
            if [ ! -f $vcf.csi ]; then
                bcftools index -c $vcf
            fi
        done
        
        # Get unique chromosome names
        CHROMOSOMES=$(awk '{print $1}' ~{query_variants} | sort | uniq)
        # Split the query_variants file per chromosome
        VCF_FILES=""
        
        for chr in $CHROMOSOMES; do
            awk -v chr=$chr '$1 == chr' ~{query_variants} > ${chr}_query_subset_variants_list.txt
            if [ -f ${chr}_query_subset_variants_list.txt ] && [ -f ~{query_samples} ]; then
                bcftools view --regions-file ${chr}_query_subset_variants_list.txt --samples-file ~{query_samples} ${chr}.dose.vcf.gz > ~{prefix}_${chr}_query_subset.vcf
            else
                bcftools view --regions-file ${chr}_query_subset_variants_list.txt ${chr}.dose.vcf.gz > ~{prefix}_${chr}_query_subset.vcf
            fi
            VCF_FILES="${VCF_FILES} ~{prefix}_${chr}_query_subset.vcf"
        done
        
        bcftools concat ${VCF_FILES} -o ~{prefix}_query_extracted.vcf
        # extract snp INFO
        python3 /scripts/extract_snp_info.py --vcf ~{prefix}_query_extracted.vcf --out ~{prefix}_query_subset_extracted_snps_info.tsv
        # extract FORMAT fields
        python3 /scripts/extract_vcf_info.py --vcf ~{prefix}_query_extracted.vcf --out ~{prefix}_extracted --extract ~{extract_item}
    >>>
    
    output {
        File vcf = prefix + "_query_extracted.vcf"
        File snp_info = prefix + "_query_subset_extracted_snps_info.tsv"
        File? gt_info = prefix + "_extracted_GT.csv"
        File? ds_info = prefix + "_extracted_DS.csv"
        File? gp_info = prefix + "_extracted_GP.csv"
    }
    
    runtime {
        docker: "docker.io/library/extract"
        memory: "32G"
        disks: "local-disk ~{disk_size_gb} HDD"
    }
}