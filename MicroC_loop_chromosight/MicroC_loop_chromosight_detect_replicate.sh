# FOR SOME REASON, THERE IS MEMORY ISSUE WHEN IT IS RAN USING SBATCH
# USED INTERACTIVE NODE
# Interactive node also doesn't work sometime depending on node
####################################################################################
conda activate chromosight

sample=("G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled")
sample=("EpiG1DMSO_pooled" "EpiG1dTAG_pooled")

resolutions=("25000" "10000" "5000")

coolDir="/athena/apostoloulab/scratch/ukl4001/data/cool_norm_pooled"
resultDir="/athena/apostoloulab/scratch/ukl4001/data/loop_chromosight"
mkdir -p $resultDir

for samp in "${sample[@]}"; do
    for res in "${resolutions[@]}"; do
        if [[ "$res" == "25000" ]]; then
            pearson_cutoff=0.3
            min_separation_value=35000
            min_dist_value=50000
            max_dist_value=2000000
            perc_undetected_value=10
        elif [[ "$res" == "10000" ]]; then
            pearson_cutoff=0.3
            min_separation_value=12000
            min_dist_value=50000
            max_dist_value=2000000
            perc_undetected_value=10
        elif [[ "$res" == "5000" ]]; then
            pearson_cutoff=0.3
            min_separation_value=6000
            min_dist_value=20000
            max_dist_value=2000000
            perc_undetected_value=10
        else
            echo "Unknown resolution: $res"
            continue
        fi
        chromosight detect "${coolDir}/${samp}_${res}bp_KR.cool" \
            "${resultDir}/${samp}_chromosight_${res}bp" \
            --threads=16 \
            --pattern=loops \
            --pearson="${pearson_cutoff}" \
            --min-separation="${min_separation_value}" \
            --min-dist="${min_dist_value}" \
            --max-dist="${max_dist_value}" \
            --perc-undetected="${perc_undetected_value}"
    done
done


####################################################################################
# Quantify chromosight loops for post-processing

sample=("G1DMSO_pooled" "G1dTAG_pooled" "G1A485_pooled" "EpiG1DMSO_pooled" "EpiG1dTAG_pooled")
resolutions=("25000" "10000" "5000")

coolDir="/athena/apostoloulab/scratch/ukl4001/data/cool_norm_pooled"
resultDir="/athena/apostoloulab/scratch/ukl4001/data/chromosight"

for res in "${resolutions[@]}"; do
    for samp in "${sample[@]}"; do
        chromosight quantify ${resultDir}/chromo_union_${res}bp.bedpe \
                        ${coolDir}/${samp}_${res}bp_KR.cool \
                        ${resultDir}/${samp}_union_${res}bp_pu100pz100 \
                        --pattern=loops \
                        --threads=64 \
                        --perc-undetected=100 \
                        --perc-zero=100
    done
done
