# convert_vcf.R
# Description: take the combined vcf file
# extract only the samples
library(SeqArray)
library(dplyr)

vcf_file <- "raw_data/malaria_gen5_biallelic_snps.vcf.gz"

vcf_header <- seqVCF_Header(vcf_file)

# recode header format for AD from R to .
vcf_header$format$Number[vcf_header$ID == "AD"] <- "."

# info columns to keep
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "ClippingRankSum",
                 "DP", "FS", "InbreedingCoeff", "MQ", "MQRankSum",
                 "QD", "ReadPosRankSum",
                 "SOR", "VQSLOD", "SNPEFF_AMINO_ACID_CHANGE",
                 "SNPEFF_CODON_CHANGE", "SNPEFF_EFFECT", "SNPEFF_EXON_ID",
                 "SNPEFF_FUNCTIONAL_CLASS", "SNPEFF_GENE_BIOTYPE", 
                 "SNPEFF_GENE_NAME", "SNPEFF_IMPACT", "SNPEFF_TRANSCRIPT_ID")

# format columns to keep
fmt.info <- c("AD", "DP", "GQ", "GT", "PL", "PGT",  "PID")



seqVCF2GDS(vcf_file, "processed_data/malaria_gen5.gds",
           header = vcf_header, info.import = info.import, fmt.import = fmt.info)
