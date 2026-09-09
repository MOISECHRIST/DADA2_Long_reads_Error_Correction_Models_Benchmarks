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
errorEstFunc.list <- c("loessErrfun"=loessErrfun, 
                      "PacBioErrfun"=PacBioErrfun, 
                      "makeBinnedQualErrfun"=makeBinnedQualErrfun,
                      "loessErrfun_mod0"=loessErrfun_mod0, 
                      "loessErrfun_mod1"=loessErrfun_mod1, 
                      "loessErrfun_mod2"=loessErrfun_mod2, 
                      "loessErrfun_mod3"=loessErrfun_mod3, 
                      "loessErrfun_mod4"=loessErrfun_mod4)
err <- list()
dd_res <- list()
track <- list()

for(func_name in names(errorEstFunc.list)){

  func <- errorEstFunc.list[[func_name]]

  message(paste0("Using the error function : ",func_name))
  err_obj <- tryCatch({
    if (func_name == "makeBinnedQualErrfun") {
      learnErrors(
        path.filts,
        errorEstimationFunction = func(c(3, 10, 17, 22, 27, 35, 40)),
        multithread = TRUE
      )
    } else {
      learnErrors(
        path.filts,
        errorEstimationFunction = func,
        multithread = TRUE
      )
    }
  }, error = function(e) {
    warning(paste0("Failed to learn errors for ", func_name, ": ", e$message))
    return(NULL)
  })
  
  if (is.null(err_obj) || is.null(err_obj$err_out)) {
    message(paste0("Skipping dada() for ", func_name, " due to error estimation failure."))
    next
  }
  
  err[[func_name]] <- err_obj

  pdf(file.path(path.figures, paste0("learn_error_plot_",func_name,".pdf")))
  print(plotErrors(err[[func_name]]))
  dev.off()

  saveRDS(err[[func_name]], file.path(path.rds, paste0("learn_error_data_",func_name,".rds")))
  dd_res[[func_name]] <- dada(path.filts, err=err[[func_name]], multithread=TRUE)

  track[[func_name]] <- sapply(dd_res[[func_name]], function(x) sum(x$denoised))
  saveRDS(dd_res[[func_name]], file.path(path.rds, paste0("dada_results_",func_name,".rds")))
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


