#!/bin/bash 

STEP=${1:-0}

DATASET_LIST=(UniBe_Revio
UniBe_Sequel
UniBe_Vega
UniL_Revio-Kinnex-Non-Segemented-SPRQ-NX
UniL_Revio-Kinnex-Segmented-SPRQ
UniL_Revio-Kinnex-Segmented-SPRQ-NX
UniL_Revio-Standard-NON-SPRQ
UniL_Revio-Standard-SPRQ
UniL_Sequel2)

NBASES_LIST=(1e+07
1e+06
1e+05
1e+04)

SEED_LIST=(1 2 3 4 5 6 7 8 9 10)

if [ $STEP -eq 0 ]; then

    for dataset in "${DATASET_LIST[@]}"; do
        sbatch scripts/bbtools_sampling.sh "data/concat/${dataset}" "1e+08" 41
    done

elif [ $STEP -eq 1 ]; then 

    for nbases in "${NBASES_LIST[@]}"; do 
        for dataset in "${DATASET_LIST[@]}"; do
            for seed in "${SEED_LIST[@]}"; do
                sbatch scripts/bbtools_sampling.sh "data_normalized/1e+08/${dataset}_41" "$nbases" "$seed"
            done 
        done
    done 

else 
    echo "Option not available"
    exit 1
fi