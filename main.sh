#!/bin/bash

##------------------------------------------------------------------------
## Author : MEKA Moise
## Email : moise.meka@students.unibe.ch
## Description : With this bash script, I will run dada2 workflow (QC, Filter, Trimming, Learn Errors, Denoising) on my dataset 
## Creation date : 08-09-2026
##------------------------------------------------------------------------

#SBATCH --partition=pibu_el8
#SBATCH --mail-user=moise.meka@students.unibe.ch
#SBATCH --mail-type=start,end,fail
#SBATCH --job-name="dada2_workflow"
#SBATCH --mem=150GB
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out

set -euo pipefail
module purge 2>/dev/null || true 

# Path to directory with all data sets
ALL_DATA_DIR=$1
DATASET_LIST=($(ls "$ALL_DATA_DIR"))
NBASES=$2
if [ -z "$NBASES" ]; then 
  NBASES="1e+08"
fi

data_dir="${DATASET_LIST[${SLURM_ARRAY_TASK_ID}]}"
echo "Working on : ${data_dir}"

apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
  Rscript ./scripts/quality_control_trimming.r "${ALL_DATA_DIR}/${data_dir}" "results/${data_dir}"

apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
  Rscript ./scripts/dada2_analysis.R "results/${data_dir}/Filtered" "results/${data_dir}" "$NBASES"
