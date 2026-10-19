#!/bin/bash

##------------------------------------------------------------------------
## Author : MEKA Moise
## Email : moise.meka@students.unibe.ch
## Description : This bash script I will run bbtools to sample fastq for 1e+08, 1e+07, 1e+06, 1e+05, 1e+04 
## Creation date : 17-10-2026
##------------------------------------------------------------------------

#SBATCH --partition=pibu_el8
#SBATCH --mail-user=moise.meka@students.unibe.ch
#SBATCH --mail-type=start,end,fail
#SBATCH --job-name="bbtools_sampling"
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out

if [ -z "$1" ]; then 
    echo "Missing parameters"
    echo "USAGE : sbatch $0 </path/to/datadir> <nbases>"
fi

DATA_DIR=$1
NBASES=$(printf "%.0f" "$2")
CONTAINER="${PWD}/containers/bbtools_40.01.sif"
OUTDIR="${PWD}/data_normalized/$2/$(basename $DATA_DIR)"

fastq_file=$(ls "$DATA_DIR")

mkdir -p "$OUTDIR"
apptainer exec --cleanenv --bind "$PWD:/workdir" \
    --pwd /workdir \
    "${CONTAINER}" \
    reformat.sh \
        in1="${DATA_DIR}/${fastq_file}" out1="${OUTDIR}/${fastq_file}" \
        samplebasestarget="${NBASES}" sampleseed=41
