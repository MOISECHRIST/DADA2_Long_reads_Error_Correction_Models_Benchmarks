#!/bin/bash
func_list=(loessErrfun PacBioErrfun makeBinnedQualErrfun loessErrfun_mod0  loessErrfun_mod1 loessErrfun_mod2 loessErrfun_mod3 loessErrfun_mod4)
RAWDATASETCOUNT=$(ls data/raw | wc -l)
step=${1:-1}
if [ $step -eq 1 ]; then
    for func_name in "${func_list[@]}"; do
        sbatch --array=0-$((RAWDATASETCOUNT-1)) main.sh data/raw "${func_name}"
    done 
elif [ $step -eq 2 ]; then
    for func_name in "${func_list[@]}"; do
        for prop_seed in $(echo "1 2 3 4 5 6 7 8 9 10 41"); do 
            DATASETCOUNT=$(ls data/"${prop_seed}" | wc -l)
            sbatch --array=0-$((DATASETCOUNT-1)) main.sh data/"${prop_seed}" "${func_name}"
        done 
    done 
else 
    for func_name in "${func_list[@]}"; do
        for nbases in $(echo "1e+04 1e+05 1e+06 1e+07"); do 
            sbatch --array=0-$((RAWDATASETCOUNT-1)) main.sh data/raw "${func_name}" "${nbases}" 
        done 
    done 
fi