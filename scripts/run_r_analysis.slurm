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
#SBATCH --job-name="prop_data_analysis"
#SBATCH --mem=150GB
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out


apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir ./containers/dada2-pipeline-1.3.sif \
    Rscript ./scripts/data_analysis_prop.R
