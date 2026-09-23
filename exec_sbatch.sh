#!/bin/bash
func_list=(loessErrfun PacBioErrfun makeBinnedQualErrfun loessErrfun_mod0  loessErrfun_mod1 loessErrfun_mod2 loessErrfun_mod3 loessErrfun_mod4)
RAWDATASETCOUNT=$(ls data/raw | wc -l)
step=${1:-0}

if [ $step -eq 0 ]; then
    #Raw data QC step
    sbatch --array=0-$((RAWDATASETCOUNT-1)) exec_quality_control.sh data/raw

    #Proportional approach QC step
    for prop_seed in $(echo "1 2 3 4 5 6 7 8 9 10 41"); do 
        DATASETCOUNT=$(ls data/"${prop_seed}" | wc -l)
        sbatch --array=0-$((DATASETCOUNT-1)) exec_quality_control.sh data/"${prop_seed}"
    done 
elif [ $step -eq 1 ]; then
    #LearnError + Denoising of the raw dataset for each function
    for func_name in "${func_list[@]}"; do
        sbatch --array=0-$((RAWDATASETCOUNT-1)) main.sh data/raw "${func_name}"
    done 
elif [ $step -eq 2 ]; then
    #LearnError + Denoising of the dataset in the proportional approach for each function
    for func_name in "${func_list[@]}"; do
        for prop_seed in $(echo "1 2 3 4 5 6 7 8 9 10 41"); do 
            DATASETCOUNT=$(ls data/"${prop_seed}" | wc -l)
            sbatch --array=0-$((DATASETCOUNT-1)) main.sh data/"${prop_seed}" "${func_name}"
        done 
    done 
elif [ $step -eq 3 ]; then
    #LearnError + Denoising of the dataset in the nbases approach for each function
    for func_name in "${func_list[@]}"; do
        for nbases in $(echo "1e+04 1e+05 1e+06 1e+07"); do 
            sbatch --array=0-$((RAWDATASETCOUNT-1)) main.sh data/raw "${func_name}" "${nbases}" 
        done 
    done 
else 
    #LearnError + Denoising of the dataset in the nbases approach for each function + randomized sampling + multiple seed
    for func_name in "${func_list[@]}"; do
        for nbases in $(echo "1e+04 1e+05 1e+06 1e+07"); do 
            for seed in $(echo "100 101 102 103"); do
                sbatch --array=0-$((RAWDATASETCOUNT-1)) main.sh data/raw "${func_name}" "${nbases}" "${seed}"
            done 
        done 
    done
fi