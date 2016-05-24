#!/bin/bash
# filter_combined.vcf
# Description: create filtered vcf file using MAF filtering on
# field isolates only. Then use that information to construct
# a GDS format file for globally polymorphic sites on the mixtures.
set -e
isolate_ids=metadata/non_field_isolates.txt
mixture_ids=metadata/moi_isolates.txt
# use vcftools to establish MAF cut-off

# vcftools --gzvcf processed_data/combined_filtered.vcf.gz \
# 	--remove $isolate_ids \
# 	--maf 0.05 \
# 	--recode \
# 	--recode-INFO-all \
# 	--stdout | bgzip -c > processed_data/global_polymorphic_snps.vcf.gz
# tabix -p vcf processed_data/global_polymorphic_snps.vcf.gz

# # extract positions
# zcat processed_data/global_polymorphic_snps.vcf.gz | grep -v ^# | cut -f1,2 > processed_data/global_snps_positions.txt

# input those plus moi_isolates and then convert to GDS format
vcftools --gzvcf processed_data/combined_filtered.vcf.gz \
	--keep $mixture_ids \
	--positions processed_data/global_snps_positions.txt \
	--recode \
	--recode-INFO-all \
	--stdout | bgzip -c > processed_data/mixtures_only.vcf.gz
tabix -p vcf processed_data/mixtures_only.vcf.gz

Rscript -e "library(SeqArray); seqVCF2GDS('processed_data/mixtures_only.vcf.gz', 'processed_data/mixtures_only.gds')"