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
SEED_LIST=(1 2 3 4 5 6 7 8 9 10 41)
for seed in "${SEED_LIST[@]}"; do 
    for dataset in "${DATASET_LIST[@]}"; do
        sbatch scripts/seqkit_sampling.slurm "data_normalized/41/${dataset}" "$seed"
    done
done 