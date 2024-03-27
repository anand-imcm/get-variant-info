# WDL Workflow for Extracting Variant Information

[![Open](https://img.shields.io/badge/Open-Dockstore-blue)](https://dockstore.org/workflows/github.com/IMCM-OX/get-variant-info:main?tab=info)&nbsp;&nbsp;
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/IMCM-OX/get-variant-info/publish.yml)&nbsp;&nbsp;
![GitHub release (with filter)](https://img.shields.io/github/v/release/IMCM-OX/get-variant-info)&nbsp;&nbsp;

> [!TIP]
> To import the workflow into your Terra workspace, click on the above Dockstore badge, and select 'Terra' from the 'Launch with' widget on the Dockstore workflow page.


This repository contains a WDL (Workflow Description Language) workflow for extracting information from a set of imputed VCF files using a list of query variants or sample IDs.

The workflow extracts the following information:

- Chromosome
- Position
- Reference allele
- Alternate allele
- Allele frequency (AF)
- Minor allele frequency (MAF)
- Imputation accuracy (R2)
- Empirical R-square (ER2)
- Genotype (GT)
- Estimated Alternate Allele Dosage (DS)
- Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 (GP)

The output is a set of files containing the extracted information.

## Workflow Inputs

- `query_variants`: A tab-delimited file with a list of query variants. Each line should be formatted as: Chromosome, Pos, ID, Ref, Alt. (required)
- `query_samples`: A file with a list of sample IDs. Each line should contain one sample ID. (optional)
- `imputed_vcf`: Array of imputed VCF files and their indices. VCF files should be in .vcf.gz format and indices in CSI or TBI format. (required)
- `prefix`: Prefix for the output files. (required)
- `extract_item`: A string specifying the information to extract from the FORMAT field of the VCF file. The available choices are GT, DS, and GP. Please provide as a comma-separated string. (required)

## Workflow Outputs

- `snp_info`: File containing SNP information.
- `genotype_info`: Optional file containing genotype information.
- `dosage_info`: Optional file containing estimated alternate allele dosage information.
- `geno_prob_info`: Optional file containing estimated posterior probabilities for genotypes.


## Components

- **Python packages**
  - pysam
  - pandas
  - argparse
- **Tools**
  - bcftools
- **Containers**
  - ghcr.io/IMCM-OX/get-variant-info
