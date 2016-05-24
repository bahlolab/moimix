#!/bin/bash
# clean_cobmine_vcfs.sh
# Description: filter vcf files for each of the nuclear chromosomes,
# combine using bcftools merge and then convert to GDS format in R.
set -e
vcfdir=./raw_data/vcf_files/
ref=Pfalciparum.genome.fasta
outdir=processed_data/
tempdir=~/tmp/
# step 1, select only SNPs + VQSLOD > 0, regiontype = core
vcf_chrom=$(find $vcfdir -regex ".*\.*[0-9][0-9]_v3.*.vcf.gz" | sort -u)

for chr in $vcf_chrom 
do
	(stem=$(basename $chr .vcf.gz)
	echo $stem
	vcftools --gzvcf $chr \
	--remove-indels \
	--min-alleles  2 \
	--max-alleles 2 \
	--remove-filtered-all \
	--recode \
	--recode-INFO-all \
	--stdout | bgzip -c > ${tempdir}${stem}.vcf.gz 
	tabix -p vcf ${tempdir}${stem}.vcf.gz) &

done
wait

filtered_vcf=$(find $tempdir -name "*.vcf.gz" | sort -u)

bcftools concat -Oz -o ${outdir}combined_filtered.vcf.gz $filtered_vcf