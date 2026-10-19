#!/bin/bash 

DATASET_LIST=(UniBe_Revio
UniBe_Sequel
UniBe_Vega
UniL_Revio-Kinnex-Non-Segemented-SPRQ-NX
UniL_Revio-Kinnex-Segmented-SPRQ
UniL_Revio-Kinnex-Segmented-SPRQ-NX
UniL_Revio-Standard-NON-SPRQ
UniL_Revio-Standard-SPRQ
UniL_Sequel2)

NBASES_LIST=(1e+08
1e+07
1e+06
1e+05
1e+04)

for nbases in "${NBASES_LIST[@]}"; do 
    for dataset in "${DATASET_LIST[@]}"; do
        sbatch scripts/bbtools_sampling.sh "data/concat/${dataset}" "$nbases"
    done
done 