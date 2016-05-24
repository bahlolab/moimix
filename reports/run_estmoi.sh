#!/bin/bash
set -e
# run estMOI with default settings to begin with
clean_bams=$(find raw_data/bam_files/mixture_isolates/ -name "*.bam" | sort -u)
fasta=Pfalciparum.genome.fasta
clean_vcf=processed_data/mixtures_only_fixed_header.vcf.gz
echo $clean_bams
echo $fasta
echo $clean_vcf
start=`date +%s`

for bam in ${clean_bams} 
do
    (stem=$(basename $bam .bam)
    if [ ! -f ${bam}.bai ]; then
    	echo "Bam index not found!"
    	samtools index ${bam}
	fi
    echo "Run estMOI on ${stem}"
    estMOI_1.03 $bam $clean_vcf $fasta --readl 100 --out=processed_data/moi_results/${stem}) &
done;
wait

end=`date +%s`
runtime=$((end-start))

echo "estMOI completed in ${runtime} seconds" 
