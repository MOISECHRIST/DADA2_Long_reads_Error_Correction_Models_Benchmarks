#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 

#Import my functions
source("scripts/functions.R")

#Check whether the number of parameters is 2
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=2){
  print("This script takes the following parameters:")
  print("   (1) Path to input filtered reads directory")
  print("   (2) Path to output directory")
  stop("Error: Required arguments are not provided.", call.=FALSE)
}

#Setup inputs and outputs 
path.input <- args[1]
path.output <- args[2]
path.rds <- file.path(path.output, "RDS")
path.figures <- file.path(path.output, "Figure")
path.filts <- list.files(path.input, pattern="fastq.gz", full.names=TRUE)
dir.create(path.figures, recursive = TRUE, showWarnings = FALSE)
dir.create(path.rds, recursive = TRUE, showWarnings = FALSE)

#Learning errors 
errorEstFunc.list <- c(loessErrfun, PacBioErrfun, makeBinnedQualErrfun,
                       loessErrfun_mod0, loessErrfun_mod1, loessErrfun_mod2, 
                       loessErrfun_mod3, loessErrfun_mod4)
err <- list()
dd_res <- list()
track <- list()
for(func in errorEstFunc.list){
  func.name <- deparse(substitute(func))
  if(func.name == "makeBinnedQualErrfun"){
    err[[func.name]] <- learnErrors(path.filts, 
                                    errorEstimationFunction=func(c(3, 10, 17, 22, 27, 35, 40)), 
                                    multithread=TRUE)
  } else {
    err[[func.name]] <- learnErrors(path.filts, errorEstimationFunction=func, 
                                    multithread=TRUE)
  }
  pdf(file.path(path.figures, paste0("learn_error_plot_",func.name,".pdf")))
  plotErrors(err[[func.name]])
  dev.off()
  saveRDS(err[[func.name]], file.path(path.rds, paste0("learn_error_data_",func.name,".rds")))
  dd_res[[func.name]] <- dada(path.filts, err=err[[func.name]], multithread=TRUE)
  track[[func.name]] <- sapply(dd_res[[func.name]], function(x) sum(x$denoised))
  saveRDS(dd_res[[func.name]], file.path(path.rds, paste0("dada_results_",func.name,".rds")))
}

track <- as.data.frame(track)
track$sample <- basename(path.filts)
track
track_long <- pivot_longer(track, cols = -sample,
                           names_to = "method", values_to = "denoised")

ggplot(track_long, aes(x = reorder(method, denoised, FUN = median), y = denoised)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, color = "darkblue") +
  labs(x = "Error model", y = "Reads denoised")

ggsave(file.path(path.figures, "track_denoised_boxplot.pdf"), width = 8, height = 6)






