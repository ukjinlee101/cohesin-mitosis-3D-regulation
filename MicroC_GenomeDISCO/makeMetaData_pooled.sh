#!/usr/bin/env bash
set -euo pipefail


BINDIR="/athena/apostoloulab/scratch/ukl4001/temp/binned.50kb"
OUTDIR="/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_GenomeDISCO/metadata"
mkdir -p $OUTDIR

samples1=(
    "G1DMSO_pooled_allRes"
    "G1dTAG_pooled_allRes"
)

samples2=(
    "G1DMSO_pooled_allRes"
    "G1A485_pooled_allRes"
)

samples3=(
    "G1DMSO_pooled_allRes"
    "EpiG1DMSO_pooled_allRes"
)

samples4=(
    "G1DMSO_pooled_allRes"
    "EpiG1dTAG_pooled_allRes"
)

samples5=(
    "EpiG1DMSO_pooled_allRes"
    "EpiG1dTAG_pooled_allRes"
)

k=58
# Loop through each pair of samples to generate the files

for samples in samples1 samples2 samples3 samples4 samples5; do
    declare -n arr="$samples"

    for (( i=0; i<${#arr[@]}; i++ )); do
        for (( j=i+1; j<${#arr[@]}; j++ )); do
            sample1=${arr[i]}
            sample2=${arr[j]}
            file1_path="${BINDIR}/${sample1}.binned.50kb.gz"
            file2_path="${BINDIR}/${sample2}.binned.50kb.gz"

            file1="${OUTDIR}/metadata_${k}.samples"
            file2="${OUTDIR}/metadata_${k}.pairs"
            # Write to file1.txt
            echo -e "${sample1}\t${file1_path}" >> "$file1"
            echo -e "${sample2}\t${file2_path}" >> "$file1"

            # Write to file2.txt
            echo -e "${sample1}\t${sample2}" >> "$file2"

            k=$((k+1))
            echo "Files generated: $file1 and $file2"
        done
    done

done
