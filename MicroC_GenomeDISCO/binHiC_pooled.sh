#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/athena/apostoloulab/scratch/ukl4001/MicroC_MtoG1_Batch123_pooled/pipeline/result/hic_pooled"
TEMPWORKDIR="/athena/apostoloulab/scratch/ukl4001/temp"
mkdir -p $TEMPWORKDIR/binned.50kb

for SAMPLE in $(find "$WORKDIR" -maxdepth 1 -type f -name "*.hic" | while read F; do basename $F | sed 's/.hic$//'; done)
do
	echo "Processing $SAMPLE"
	TEMP_FILE="$TEMPWORKDIR/temp.txt"
	CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX")
	for chr in "${CHROMOSOMES[@]}"
	do
		echo "Processing $chr"
		java -Xmx256g -Djava.awt.headless=true -jar /athena/apostoloulab/scratch/ukl4001/001_software/juicer_tools_1.22.01.jar dump \
					observed NONE \
					${WORKDIR}/${SAMPLE}.hic \
					${chr} ${chr} BP 50000 \
					${TEMP_FILE}
		awk -v chr=$chr '{print chr, $1, chr, $2, $3}' OFS='\t' $TEMP_FILE >> "${TEMPWORKDIR}/binned.50kb/${SAMPLE}.binned.50kb"

	done
	gzip "${TEMPWORKDIR}/binned.50kb/${SAMPLE}.binned.50kb"
	rm $TEMP_FILE
done

WORKDIR="/athena/apostoloulab/scratch/ukl4001/MicroC_MtoG1_Batch4_pooled/pipeline/result/hic_pooled"
TEMPWORKDIR="/athena/apostoloulab/scratch/ukl4001/temp"
mkdir -p $TEMPWORKDIR/binned.50kb

for SAMPLE in $(find "$WORKDIR" -maxdepth 1 -type f -name "*.hic" | while read F; do basename $F | sed 's/.hic$//'; done)
do
	echo "Processing $SAMPLE"
	TEMP_FILE="$TEMPWORKDIR/temp.txt"
	CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chrX")
	for chr in "${CHROMOSOMES[@]}"
	do
		echo "Processing $chr"
		java -Xmx64g -Djava.awt.headless=true -jar /athena/apostoloulab/scratch/ukl4001/001_software/juicer_tools_1.22.01.jar dump \
					observed NONE \
					${WORKDIR}/${SAMPLE}.hic \
					${chr} ${chr} BP 50000 \
					${TEMP_FILE}
		awk -v chr=$chr '{print chr, $1, chr, $2, $3}' OFS='\t' $TEMP_FILE >> "${TEMPWORKDIR}/binned.50kb/${SAMPLE}.binned.50kb"

	done
	gzip "${TEMPWORKDIR}/binned.50kb/${SAMPLE}.binned.50kb"
	rm $TEMP_FILE
done
