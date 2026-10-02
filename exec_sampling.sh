#!/bin/bash 

DATASET_LIST=(LIB_16S_KINNEX-NON-SEGMENTED_REVIO_SPRQ_NX  LIB_16S_KINNEX_SEGMENTED_REVIO_SPRQ_NX  LIB_16S_STANDARD_REVIO_SPRQ  Revio_UniBe  LIB_16S_KINNEX_SEGMENTED_REVIO_SPRQ  LIB_16S_STANDARD_REVIO_NON_SPRQ  LIB_16S_STANDARD_Sequel2_NON_SPRQ  Sequel_UniBe)

for seed in "1 2 3 4 5 6 7 8 9 10 41"; do 
    for dataset in "${DATASET_LIST[@]}"; do
        sbatch scripts/seqkit_sampling.slurm "data_normalized/41/${dataset}" "$seed"
    done
done 