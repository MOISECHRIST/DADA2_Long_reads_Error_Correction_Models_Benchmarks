#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr") 
library(phyloseq, verbose = FALSE, quietly = TRUE); packageVersion("phyloseq")
library(Biostrings, verbose = FALSE, quietly = TRUE); packageVersion("Biostrings")

#Define functions
readMultiRDS <- function(list_of_paths){
  res <- list()
  for(file_path in list_of_paths){
    ref_name <- paste0(basename(dirname(dirname(file_path))), "_", basename(file_path))
    res[[ref_name]] <- readRDS(file_path)
  }
  return(res)
}

computeDistances <- function(dada2_results_data, dataset_list, errFunc_list, dataset_name){
  distances <- list()
  for(item1 in errFunc_list){
    tmp <- c()
    ref_names <- c()
    full <- getErrors(dada2_results_data[[paste0(dataset_name,"_dada_results_",item1)]])
    for(item2 in dataset_list){
      if(endsWith(item2, item1)){
        ref_names <- c(ref_names, 
                       sapply(strsplit(item2 , "_dada_results_"), function(x) x[1]))
        tmp <- c(tmp, 
                 wavethresh::l2norm(full, getErrors(dada2_results_data[[item2]])))
      }
    }
    tmp <- as.matrix(tmp)
    rownames(tmp) <- ref_names
    distances[[item1]] <- tmp
  }
  distances <- as.data.frame(distances)
  dataset.prop <- sapply(strsplit(rownames(distances),"_"), function(x) x[3])
  dataset.prop[is.na(dataset.prop)] <- "100"
  distances$dataset.prop <- as.numeric(dataset.prop)
  
  used.seed <- sapply(strsplit(rownames(distances),"_"), function(x) x[4])
  distances$used.seed <- as.numeric(used.seed)
  return(distances)
}

transToLong.prop <- function(distances){
  distances.long <- distances |>
    pivot_longer(cols = ends_with(".rds"),
                 values_to = "distances",
                 names_to = "errors.function") |>
    mutate(
      errors.function = sub(".rds","",errors.function),
      dataset.prop = factor(dataset.prop, levels = c("1", "5", "10", "25", "50", "100"))
    )
  return(distances.long)
}

plotDistancePoint.prop <- function(distances.long){
  distances.long |>
    ggplot()+
    geom_point(aes(x=dataset.prop, y=distances), alpha=0.4)+
    facet_grid(cols=vars(errors.function))+
    labs(x="Dataset proportion (%)",
         y="Distance to the full dataset") + theme_bw()
}

plotDistanceBoxplot.prop <- function(distances.long){
  distances.long |>
    ggplot()+
    geom_boxplot(aes(x=dataset.prop, y=distances), alpha=0.4)+
    facet_grid(cols=vars(errors.function))+
    labs(x="Dataset proportion (%)",
         y="Distance to the full dataset") + theme_bw()
}

#Input data
path.alldataset <- file.path(list.files("results", full.names = T), "RDS")
paths_learn_errors_data <- list.files(path.alldataset, pattern = "learn_error", full.names = T) 
paths_dada2_results_data <- list.files(path.alldataset, pattern = "dada_results", full.names = T) 
ref.db <- file.path("refSeq/SILVA-v138.2-16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz")

results.path <- "summary"
dir.create(results.path)

#Load data
learn_errors_data <- readMultiRDS(paths_learn_errors_data)
dada2_results_data <- readMultiRDS(paths_dada2_results_data)

#Sequel_UniBe Distances
sequel_unibe_dataset <- names(dada2_results_data)[startsWith(names(dada2_results_data),"Sequel_UniBe")]
sequel_unibe_distances <- computeDistances(dada2_results_data, 
                                     sequel_unibe_dataset, 
                                     c("loessErrfun_mod0.rds", 
                                        "PacBioErrfun.rds", 
                                       "loessErrfun.rds"),
                                     "Sequel_UniBe")

sequel_unibe_distances.long <- transToLong.prop(sequel_unibe_distances)
plotDistancePoint.prop(sequel_unibe_distances.long)
ggsave(file.path(results.path,"Sequel_UniBe_dataset_distance_plot.pdf"))
plotDistanceBoxplot.prop(sequel_unibe_distances.long)
ggsave(file.path(results.path,"Sequel_UniBe_dataset_distance_boxplot.pdf"))

#Revio_UniBe Distances 
revio_unibe_dataset <- names(dada2_results_data)[startsWith(names(dada2_results_data),"Revio_UniBe")]
revio_unibe_distances <- computeDistances(dada2_results_data, 
                                     revio_unibe_dataset, 
                                     c("loessErrfun.rds", 
                                       "PacBioErrfun.rds", 
                                       "makeBinnedQualErrfun.rds",
                                       "loessErrfun_mod0.rds"),
                                     "Revio_UniBe")

revio_unibe_distances.long <- transToLong.prop(revio_unibe_distances)
plotDistancePoint.prop(revio_unibe_distances.long)
ggsave(file.path(results.path,"Revio_UniBe_dataset_distance_plot.pdf"))
plotDistanceBoxplot.prop(revio_unibe_distances.long)
ggsave(file.path(results.path,"Revio_UniBe_dataset_distance_boxplot.pdf"))

#Sequence table
seq_tables <- list()
taxaAssign <- list()
N <- length(names(dada2_results_data))
n=1
for (ref_name in names(dada2_results_data)){
  cat("[",n,"/",N,"] : ", ref_name,"\n")
  seq_tables[[ref_name]] <- makeSequenceTable(dada2_results_data[[ref_name]])
  n<-n+1
}
saveRDS(seq_tables,file.path(results.path, "all_dataset_sequence_table.rds"))

#Assign taxonomy
n=1
for (ref_name in names(dada2_results_data)){
  cat("[",n,"/",N,"] : ", ref_name,"\n")
  seqtab.nochim <- removeBimeraDenovo(seq_tables[[ref_name]], method="consensus", multithread=TRUE, verbose=TRUE)
  taxaAssign[[ref_name]] <- assignTaxonomy(seqtab.nochim, ref.db,multithread = T)
  n<-n+1
}
saveRDS(taxaAssign,file.path(results.path, "all_dataset_taxonomy_assignment.rds"))





