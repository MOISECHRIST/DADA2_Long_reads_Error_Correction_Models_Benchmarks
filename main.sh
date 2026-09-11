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
#SBATCH --time=72:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out

set -euo pipefail
module purge 2>/dev/null || true 

# Path to directory with all data sets
ALL_DATA_DIR=$1
if [-z "$ALL_DATA_DIR" ]; then
  echo "ERROR : Missing parameter"
  echo "USAGE : sbatch $0 /path/to/all/dataset [nbase] [seed]"
  echo "/path/to/all/dataset : A directory with all datasets \n(one directory per dataset eg. /path/to/data where $(ls /path/to/data) gives Revio_UniBe Sequel_UniBe)."
  echo "nbase : (Optional) Number bases to consider for learning error step (DEFAULT 1e+08)."
  echo "seed : (Optional) Seed for random choice." 
  exit 1
fi
DATASET_LIST=($(ls "$ALL_DATA_DIR"))
NBASES=$2
SEED=$3

data_dir="${DATASET_LIST[${SLURM_ARRAY_TASK_ID}]}"
echo "Working on : ${data_dir}"

apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
  Rscript ./scripts/quality_control_trimming.r "${ALL_DATA_DIR}/${data_dir}" "results/${data_dir}"

if [ -z "$NBASES" ]; then 
  apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
  Rscript ./scripts/dada2_analysis.R "results/${data_dir}/Filtered" "results/${data_dir}"
elif [ -z "$SEED" ]; then 
  apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
    Rscript ./scripts/dada2_analysis.R "results/${data_dir}/Filtered" "results/${data_dir}_${NBASES}_${SEED}" "$NBASES" "$SEED"
else 
  apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
    Rscript ./scripts/dada2_analysis.R "results/${data_dir}/Filtered" "results/${data_dir}_${NBASES}" "$NBASES"
fi
