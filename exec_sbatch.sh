#!/bin/bash

for func_name in $(echo "loessErrfun PacBioErrfun makeBinnedQualErrfun loessErrfun_mod0  loessErrfun_mod1 loessErrfun_mod2 loessErrfun_mod3 loessErrfun_mod4"); do
    DATASETCOUNT=$(ls data/raw | wc -l)
    sbatch --array=0-$((DATASETCOUNT-1)) main.sh data/raw "${func_name}"
done 

for func_name in $(echo "loessErrfun PacBioErrfun makeBinnedQualErrfun loessErrfun_mod0  loessErrfun_mod1 loessErrfun_mod2 loessErrfun_mod3 loessErrfun_mod4"); do
    for prop_seed in $(echo "1 2 3 4 5 6 7 8 9 10 41"); do 
        DATASETCOUNT=$(ls data/"${prop_seed}" | wc -l)
        sbatch --array=0-$((DATASETCOUNT-1)) main.sh data/raw "${func_name}"
    done 
done 

for func_name in $(echo "loessErrfun PacBioErrfun makeBinnedQualErrfun loessErrfun_mod0  loessErrfun_mod1 loessErrfun_mod2 loessErrfun_mod3 loessErrfun_mod4"); do
    for nbases in $(echo "1e+04 1e+05 1e+06 1e+07"); do 
        N=$(ls data/raw | wc -l)
        sbatch --array=0-$((N-1)) main.sh data/raw "${func_name}" "${nbases}" 
    done 
done 