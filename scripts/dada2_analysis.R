#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 

#Import my functions
source("scripts/functions.R")

#Check whether the number of parameters
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

# Read a field (in kB) from /proc/self/status, e.g. VmHWM (peak RSS), VmRSS (current)
read_status_kb <- function(field) {
  l <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
  x <- grep(paste0("^", field, ":"), l, value = TRUE)
  if (length(x) == 0) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", x))
}

# Reset the peak-RSS counter so each step gets its own peak
reset_peak_rss <- function() {
  try(writeLines("5", "/proc/self/clear_refs"), silent = TRUE)
}

resource_log <- list()

# Run an expression and record time, CPU and RAM usage
track_resources <- function(name, expr) {
  reset_peak_rss()
  t0 <- Sys.time(); p0 <- proc.time()
  res <- force(expr)
  t1 <- Sys.time(); p1 <- proc.time()

  wall    <- as.numeric(difftime(t1, t0, units = "secs"))
  cpu_sec <- sum((p1 - p0)[c("user.self", "sys.self")])

  n_alloc <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
  resource_log[[name]] <<- data.frame(
    process_name      = name,
    start_time        = t0,
    end_time          = t1,
    wall_sec          = wall,
    cpu_sec           = cpu_sec,
    avg_cores_used    = cpu_sec / wall,                         # mean number of busy cores
    threads_available = RcppParallel::defaultNumThreads(),      # threads dada2 can use
    cpus_allocated    = n_alloc,                                # NA if not under SLURM
    peak_ram_gb       = read_status_kb("VmHWM") / 1024^2,
    end_ram_gb        = read_status_kb("VmRSS") / 1024^2
  )
  res
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

err_obj <- track_resources("learnError", tryCatch({
  n_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  if (is.na(n_cpus)) n_cpus <- TRUE
  if (func_name == "makeBinnedQualErrfun") {
    learnErrors(path.filts,
                errorEstimationFunction = func(c(3, 10, 17, 22, 27, 35, 40)),
                nbases = nbases, multithread = n_cpus, randomize = randomize)
  } else {
    learnErrors(path.filts, errorEstimationFunction = func,
                nbases = nbases, multithread = n_cpus, randomize = randomize)
  }
}, error = function(e) {
  warning(paste0("Failed to learn errors for ", func_name, ": ", e$message))
  NULL
}))

 if (is.null(err_obj) || is.null(err_obj$err_out)) {
     write.csv(do.call(rbind, resource_log),
               file.path(path.rds, paste0("execution_time_", func_name, ".csv")),
               row.names = FALSE)
     stop(paste0("Error estimation failed for ", func_name))
   }

saveRDS(err_obj, file.path(path.rds, paste0("learn_error_data_", func_name, ".rds")))

dd_res <- track_resources("Denoising",
                          dada(path.filts, err = err_obj, multithread = TRUE))
saveRDS(dd_res, file.path(path.rds, paste0("dada_results_", func_name, ".rds")))

res_df <- do.call(rbind, resource_log)
write.csv(res_df, file.path(path.rds, paste0("execution_time_", func_name, ".csv")),
          row.names = FALSE)

message("Done.")
quit(save = "no", status = 0)