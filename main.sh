#!/bin/bash

WORKING_DIR=$PWD
for data_dir in $(ls "${WORKING_DIR}/data"); do 
  echo "Working on : data/$data_dir"
  Rscript ./scripts/quality_control_trimming.r "data/$data_dir" "results/$data_dir"
  Rscript ./scripts/dada2_analysis.R "results/${data_dir}/Filtered" "results/$data_dir"
done
