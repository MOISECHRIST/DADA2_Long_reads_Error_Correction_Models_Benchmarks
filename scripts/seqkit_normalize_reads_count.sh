#!/bin/bash

##------------------------------------------------------------------------
## Author : MEKA Moise
## Email : moise.meka@students.unibe.ch
## Description : This bash script I will run seqkit to sample reads to 50%, 25%, 10%, 5% and 1%
## Creation date : 08-09-2026
##------------------------------------------------------------------------

#SBATCH --partition=pibu_el8
#SBATCH --mail-user=moise.meka@students.unibe.ch
#SBATCH --mail-type=start,end,fail
#SBATCH --job-name="seqkit_sampling"
#SBATCH --mem=150GB
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out

set -euo pipefail

# Load required module
module load SeqKit/2.6.1

# Define variables and constants
DATA_DIR=$1
SEED=${2:-41}
THREADS=$SLURM_CPUS_PER_TASK
OUTDIR="data_normalized/${SEED}"
MIN_COUNT=$(seqkit stats -j "$THREADS" "${DATA_DIR}/*/*" | awk 'NR>1 {print $4}' | sed 's/,//g' | sort -n | head -n 1)

for fastq_file in "$DATA_DIR"/*/*.fastq*; do
    dataset=$(basename $(dirname $fastq_file)) 
    mkdir -p "${OUTDIR}/${dataset}"
    out_file="${OUTDIR}/${dataset}/$(basename "$fastq_file")"
    seqkit sample --threads "$THREADS" --rand-seed "$SEED" --number "$MIN_COUNT" "$fastq_file" | gzip > "$out_file"
done