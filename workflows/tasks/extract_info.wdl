version 1.0

task extract {

    input {
        File query_variants
        File? query_samples
        Array [File] imputed_vcf
        String prefix
        String extract_item
        String docker
    }
    
    Int disk_size_gb = ceil(size(imputed_vcf, "GiB")) + 5
    
    command <<<
        for vcf in ~{sep=' ' imputed_vcf}; do
            ln -s $vcf $(basename $vcf)
            # Check if the file has the .vcf.gz extension
            if [[ $vcf == *.vcf.gz ]]; then
                # If the index file does not exist, create a csi index
                if [ ! -f $vcf.csi ]; then
                    bcftools index -c $vcf
                fi
            fi
        done
        
        # Get unique chromosome names from the query_variants file
        CHROMOSOMES=$(awk '{print $1}' ~{query_variants} | sort | uniq)
        
        # Split the query_variants file per chromosome and create a subset vcf
        VCF_FILES=""
        for chr in $CHROMOSOMES; do
            awk -v chr=$chr '$1 == chr' ~{query_variants} > ${chr}_query_subset_variants_list.txt
            if [ -f ${chr}_query_subset_variants_list.txt ]; then
                if [ -n "~{query_samples}" ] && [ -f ~{query_samples} ]; then
                    # If the list of variants and the list of samples exist, create a subset VCF
                    bcftools view --regions-file ${chr}_query_subset_variants_list.txt --samples-file ~{query_samples} ${chr}.dose.vcf.gz > ~{prefix}_${chr}_query_subset.vcf
                else
                    # If the list of samples does not exist, create a subset VCF without filtering samples
                    bcftools view --regions-file ${chr}_query_subset_variants_list.txt ${chr}.dose.vcf.gz > ~{prefix}_${chr}_query_subset.vcf
                fi
            fi
            VCF_FILES="${VCF_FILES} ~{prefix}_${chr}_query_subset.vcf"
        done
        
        # Combine the chromosome-wise subset VCFs into a single VCF
        bcftools concat ${VCF_FILES} -o ~{prefix}_query_extracted.vcf
        
        plink2 --vcf ~{prefix}_query_extracted.vcf --recode compound-genotypes --out ~{prefix}_query_extracted
        
        # extract snp INFO and FORMAT fields from the VCF
        python3 /scripts/extract_vcf_info.py --vcf ~{prefix}_query_extracted.vcf --out ~{prefix}_extracted --extract ~{extract_item}
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