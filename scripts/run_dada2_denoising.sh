#!/bin/bash

##------------------------------------------------------------------------
## Author : MEKA Moise
## Description : Cleaning + denoising (+ sequence tables) of one dataset with every 
##               learned error of that dataset (error function x nbases).
##               One array task = one dataset of READSDIR.
##------------------------------------------------------------------------

#SBATCH --partition=pshort_el8
#SBATCH --mail-user=moise.meka@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --job-name="dada2_denoising_taxonomy"
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=20
#SBATCH --time=02:00:00
#SBATCH --error=.log/errors/%x_%j.err
#SBATCH --output=.log/output/%x_%j.out

set -euo pipefail

step=${1:-0}

if [ $step -eq 0 ]; then 
  LEARNERRORS=results_tmp
  READSDIR=raw
  OUTDIR=taxa_out_nbase

  apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
    Rscript scripts/dada2_denoising-taxonomy_nbase.R $LEARNERRORS $READSDIR $OUTDIR

elif [ $step -eq 1 ]; then 
  RUNDIR=results
  OUTDIR=taxa_out_prop

  apptainer exec --cleanenv --bind "$PWD:/workdir" --pwd /workdir containers/dada2-pipeline.sif \
  Rscript scripts/dada2_denoising-taxonomy_prop.R $RUNDIR $OUTDIR
else
  echo "Option not available"
  exit 1
fi