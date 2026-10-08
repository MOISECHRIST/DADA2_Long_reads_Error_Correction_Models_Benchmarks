#!/bin/bash

##------------------------------------------------------------------------
## Author : MEKA Moise
## Email : moise.meka@students.unibe.ch
## Description : This bash script I will run seqkit to normalize datasets to get the same amount of reads per fastq file
## Creation date : 08-09-2026
##------------------------------------------------------------------------

#SBATCH --partition=pshort_el8
#SBATCH --mail-user=moise.meka@students.unibe.ch
#SBATCH --mail-type=start,end,fail
#SBATCH --job-name="seqkit_normalize"
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out

set -euo pipefail

# Load required module
module load SeqKit/2.6.1

# Define variables and constants
DATA_DIR=$1
SEED=${2:-41}
THREADS=$SLURM_CPUS_PER_TASK
OUTDIR="/data/users/${USER}/research_project/data_normalized/${SEED}"
mkdir -p "$OUTDIR"

FILES=("${DATA_DIR}"/*/*.fastq*)
if [ ! -f "${OUTDIR}/stats_before.tsv" ]; then 
    seqkit stats -T -j "$THREADS" "${FILES[@]}" > "${OUTDIR}/stats_before.tsv"
    MIN_COUNT=$(awk 'NR>1 {print $4}' "${OUTDIR}/stats_before.tsv" | sed 's/,//g' | sort -n | head -n 1)
    echo "Files: ${#FILES[@]} | Min reads: ${MIN_COUNT} | Seed: ${SEED}"
else 
    MIN_COUNT=7921
fi

for fastq_file in "${FILES[@]}"; do
    dataset=$(basename "$(dirname "$fastq_file")")
    base=$(basename "$fastq_file")
    base=${base%.gz}; base=${base%.fastq}
    mkdir -p "${OUTDIR}/${dataset}"
    seqkit sample -2 -j "$THREADS" -s "$SEED" -n "$MIN_COUNT" \
        "$fastq_file" -o "${OUTDIR}/${dataset}/${base}.fastq.gz"
done

seqkit stats -T -j "$THREADS" "${OUTDIR}"/*/*.fastq.gz > "${OUTDIR}/stats_after.tsv"