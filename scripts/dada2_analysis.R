#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 

#Import my functions
source("scripts/functions.R")

#Check whether the number of parameters is 2
args <- commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  print("This script takes the following parameters:")
  print("   (1) Path to input filtered reads directory")
  print("   (2) Path to output directory")
  print("   (3) The name of the estimation error function. (DEFAULT loessErrfun)")
  print("   (4) Number bases to consider for learning error step (DEFAULT 1e+08).")
  print("   (5) Seed for random choice.")
  stop("Error: Required arguments are not provided.", call.=FALSE)
}

#Setup inputs and outputs 
path.input <- args[1]
path.output <- args[2]
func_name <- args[3]
if (is.na(func_name)){
  func_name <- "loessErrfun"
}

nbases <- as.numeric(args[4])
if (is.na(nbases)) {
  nbases <- 1e+08
}
seed<-as.numeric(args[5])
if (is.na(seed)) {
  randomize=FALSE
} else {
  randomize=TRUE
}
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


func <- errorEstFunc.list[[func_name]]

message(paste0("Using the error function : ",func_name))
if (! is.na(seed)){
  set.seed(seed)
}
err_obj <- tryCatch({
  if (func_name == "makeBinnedQualErrfun") {
    learnErrors(
      path.filts,
      errorEstimationFunction = func(c(3, 10, 17, 22, 27, 35, 40)),
      nbases = nbases,
      multithread = TRUE,
      randomize=randomize
    )
  } else {
    learnErrors(
      path.filts,
      errorEstimationFunction = func,
      nbases = nbases,
      multithread = TRUE, 
      randomize=randomize
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

saveRDS(err_obj, file.path(path.rds, paste0("learn_error_data_",func_name,".rds")))
dd_res <- dada(path.filts, err=err_obj, multithread=TRUE)

saveRDS(dd_res, file.path(path.rds, paste0("dada_results_",func_name,".rds")))

message("Done.")
quit(save = "no", status = 0)