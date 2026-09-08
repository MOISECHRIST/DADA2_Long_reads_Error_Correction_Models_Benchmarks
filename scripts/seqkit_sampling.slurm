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
#SBATCH --array=0-4

set -euo pipefail

# Load required module
module load SeqKit/2.6.1

# Define variables and constants
DATA_DIR="${1:-}"
prop_list=(50 25 10 5 1)

if [ -z "$DATA_DIR" ]; then
    echo "ERROR : Missing parameter"
    echo "USAGE : sbatch $0 path/to/data/directory"
    echo "You have to provide the directory on which there are fastq files"
    exit 1
fi

pct="${prop_list[$SLURM_ARRAY_TASK_ID]}"
prop=$(echo "scale=2; ${pct}/100" | bc -l)
prop="0${prop}"  

OUT_DIR="$(dirname "$DATA_DIR")/$(basename "$DATA_DIR")_${pct}"
mkdir -p "$OUT_DIR"

for fastq_file in "$DATA_DIR"/*.fastq*; do
    out_file="${OUT_DIR}/$(basename "$fastq_file")"
    seqkit sample --rand-seed 41 --proportion "$prop" "$fastq_file" | gzip > "$out_file"
done