## Project structure

```         
├── Dockerfile
├── exec_quality_control.sh
├── exec_sampling.sh
├── exec_sbatch_normalized.sh
├── exec_sbatch.sh
├── main.sh
├── research_project.Rproj
├── README.md
├── refSeq/
├── results/ -> Job results on IBU cluster 
├── results_nbases/ -> Results from nbases approach
├── results_prop/ -> Results from proportion approach
├── scripts/
|   ├── dada2_analysis.R
|   ├── data_analysis_nbases.R
|   ├── data_analysis_prop.R
|   ├── functions.R
|   ├── gather_logs.py
|   ├── phyloseq_visualization.R
|   ├── quality_control_trimming.r
|   ├── run_r_analysis.slurm
|   ├── seqkit_normalize_reads_count.sh
|   └── seqkit_sampling.slurm
└── summary/
    ├── Dataset_Summary-Before_Normalization.csv
    ├── nbases/
    ├── proportion/
    ├── results_summary_learning.csv
    ├── taxonomic_distrib_top_100_all_datasets_plot.pdf
    └── taxonomic_distrib_top_100_nbases_all_datasets_plot.pdf
```


