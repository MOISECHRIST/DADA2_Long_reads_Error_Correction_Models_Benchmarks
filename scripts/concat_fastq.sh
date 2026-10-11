#!/bin/bash

##------------------------------------------------------------------------
## Author : MEKA Moise
## Email : moise.meka@students.unibe.ch
## Description : 
## Creation date : 11-10-2026
##------------------------------------------------------------------------

#SBATCH --partition=pshort_el8
#SBATCH --mail-user=moise.meka@students.unibe.ch
#SBATCH --mail-type=start,end,fail
#SBATCH --job-name="concat_fastq"
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --error=/data/users/%u/research_project/.log/errors/%x_%j.err
#SBATCH --output=/data/users/%u/research_project/.log/output/%x_%j.out

DATA_DIR=$1
OUT_PATH=$2

zcat "${DATA_DIR}/*.fastq.gz" | gzip > "$OUT_PATH"