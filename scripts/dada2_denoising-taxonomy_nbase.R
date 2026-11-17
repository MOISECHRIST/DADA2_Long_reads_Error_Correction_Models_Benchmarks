#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr")

#Check the number of parameters
args <- commandArgs(trailingOnly=TRUE)

if(length(args)<3){
  print("This script takes the following parameters:")
  print("   (1) Path to the directory with the learned errors (<run>/RDS/learn_error_data_<func>.rds)")
  print("   (2) Path to the directory with the raw reads (one directory per dataset, fastq.gz)")
  print("   (3) Path to output directory")
  print("   (4) Path to the taxonomy reference database (DEFAULT SILVA v138.2 trainset)")
  print("   (5) Path to the directory of cleaned reads (DEFAULT filtered_reads)")
  stop("Error: Required arguments are not provided.", call.=FALSE)
}

path.errors   <- args[1]
path.fastq    <- args[2]
path.output   <- args[3]
ref.db        <- if (length(args) >= 4) args[4] else
  "refSeq/SILVA-v138.2-16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
path.filtered <- if (length(args) >= 5) args[5] else "filtered_reads"

if(!dir.exists(path.errors)) stop("Directory not found: ", path.errors)
if(!dir.exists(path.fastq))  stop("Directory not found: ", path.fastq)
if(!file.exists(ref.db))     stop("Reference database not found: ", ref.db)
dir.create(path.output, recursive = TRUE, showWarnings = FALSE)

n_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cpus)) n_cpus <- TRUE

#Constants
prop.levels <- c("1e+04", "1e+05", "1e+06", "1e+07", "1e+08")
full.prop <- "1e+08"   # number of bases of the runs without "_<nbases>" suffix (default)
## Run directory: <dataset>[_<nbases>[_<seed>]]
nbases.pat <- "^(.+?)_(1e\\+0[0-9])(?:_([0-9]+))?$"


##------------------------------------------------------------------------
## 1. Dataset table : dataset x error function x nbases (x fastq file)
##------------------------------------------------------------------------
err.paths <- list.files(path.errors, pattern = "^learn_error_data_.*\\.rds$",
                        recursive = TRUE, full.names = TRUE)
if(length(err.paths) == 0) stop("No learn_error_data_*.rds found in ", path.errors)

run        <- basename(dirname(dirname(err.paths)))
has.suffix <- grepl(nbases.pat, run, perl = TRUE)
seed <- ifelse(has.suffix, sub(nbases.pat, "\\3", run, perl = TRUE), NA_character_)
seed[!is.na(seed) & seed == ""] <- NA_character_

runs.table <- data.frame(
  dataset          = ifelse(has.suffix, sub(nbases.pat, "\\1", run, perl = TRUE), run),
  error_function   = sub("^learn_error_data_", "", sub("\\.rds$", "", basename(err.paths))),
  nbases           = ifelse(has.suffix, sub(nbases.pat, "\\2", run, perl = TRUE), full.prop),
  seed             = seed,
  run              = run,
  learn_error_rds  = normalizePath(err.paths),
  stringsAsFactors = FALSE
)

## Raw reads of a run: <fastq root>/<dataset>_<nbases> if it exists (reads subsampled per nbases), 
## otherwise <fastq root>/<dataset> (shared by all the nbases of the dataset)
findFastqDir <- function(run, dataset){
  for(cand in file.path(path.fastq, c(run, dataset))){
    if(dir.exists(cand)) return(normalizePath(cand))
  }
  NA_character_
}
runs.table$fastq_dir <- vapply(seq_len(nrow(runs.table)),
                               function(i) findFastqDir(runs.table$run[i], runs.table$dataset[i]),
                               character(1))
if(anyNA(runs.table$fastq_dir)){
  message("[skipped] no raw reads directory found for: ",
          paste(unique(runs.table$run[is.na(runs.table$fastq_dir)]), collapse = ", "))
  runs.table <- runs.table[!is.na(runs.table$fastq_dir), ]
}

## One row per fastq file
datasets.table <- bind_rows(lapply(seq_len(nrow(runs.table)), function(i){
  files <- list.files(runs.table$fastq_dir[i], pattern = "\\.fastq\\.gz$", full.names = TRUE)
  if(length(files) == 0){
    message("[skipped] no fastq.gz in ", runs.table$fastq_dir[i])
    return(NULL)
  }
  rows <- runs.table[rep(i, length(files)), ]
  rows$fastq <- normalizePath(files)
  rows
}))
if(nrow(datasets.table) == 0) stop("No fastq.gz file found in ", path.fastq)

## Where the cleaned reads go: filtered_reads/<dataset>[_<nbases>]/Filtered/
## (the same raw reads are cleaned once, whatever the number of bases of the run)
datasets.table$filtered_dir   <- file.path(path.filtered, basename(datasets.table$fastq_dir))
datasets.table$filtered_fastq <- file.path(datasets.table$filtered_dir, "Filtered",
                                           basename(datasets.table$fastq))
datasets.table <- datasets.table |> arrange(dataset, nbases, error_function, fastq)
rownames(datasets.table) <- NULL
write.csv(datasets.table, file.path(path.output, "datasets_table.csv"), row.names = FALSE)
message(n_distinct(datasets.table$dataset), " datasets, ",
        n_distinct(datasets.table$error_function), " error functions, ",
        n_distinct(datasets.table$nbases), " nbases: ",
        n_distinct(datasets.table[c("run", "error_function")]), " runs to denoise")


##------------------------------------------------------------------------
## 2. Cleaning of the fastq files (skipped when the cleaned file already exists)
##------------------------------------------------------------------------
to.clean <- distinct(datasets.table, fastq, filtered_dir)
for(path.out in unique(to.clean$filtered_dir)){
  fastq.input <- to.clean$fastq[to.clean$filtered_dir == path.out]
  filts <- file.path(path.out, "Filtered", basename(fastq.input))
  todo  <- !file.exists(filts)
  if(!any(todo)){
    message("[skip] already cleaned: ", path.out)
    next
  }
  message("[cleaning] ", path.out, " (", sum(todo), " file(s))")
  dir.create(file.path(path.out, "Filtered"), recursive = TRUE, showWarnings = FALSE)
  
  #Cleaned in a temporary directory then moved: an interrupted job never leaves a 
  #truncated file that would be taken for a finished one
  path.tmp <- file.path(path.out, "Filtered_tmp")
  unlink(path.tmp, recursive = TRUE)
  dir.create(path.tmp, recursive = TRUE, showWarnings = FALSE)
  tmp.filts <- file.path(path.tmp, basename(fastq.input[todo]))
  
  track <- filterAndTrim(fastq.input[todo], tmp.filts, minQ=3, minLen=1000, maxLen=1600,
                         maxN=0, rm.phix=FALSE, maxEE=2, multithread = n_cpus)
  track <- as.data.frame(track)
  track$sample <- rownames(track)
  rownames(track) <- NULL 
  
  ok <- file.exists(tmp.filts)   # no file when no read passes the filter
  file.rename(tmp.filts[ok], filts[todo][ok])
  unlink(path.tmp, recursive = TRUE)
  
  track.file <- file.path(path.out, "track_filter-trim.csv")
  if(file.exists(track.file)){
    old <- read.csv(track.file)
    track <- bind_rows(old[!old$sample %in% track$sample, ], track)
  }
  write.csv(track, track.file, row.names = FALSE)
}

missing <- !file.exists(datasets.table$filtered_fastq)
if(any(missing)){
  message("[skipped] no cleaned reads (no read passed the filter) for: ",
          paste(unique(datasets.table$run[missing]), collapse = ", "))
  datasets.table <- datasets.table[!missing, ]
}
if(nrow(datasets.table) == 0) stop("No cleaned reads available")


##------------------------------------------------------------------------
## 3. Denoising : list of results, one element per dataset x nbases x error function
##------------------------------------------------------------------------
runs <- split(datasets.table, paste0(datasets.table$run, "_", datasets.table$error_function))
dada_results <- list()
failed <- list()
N <- length(runs)
n <- 1
for(ref_name in names(runs)){
  cat("[", n, "/", N, "] denoising : ", ref_name, "\n")
  n <- n + 1
  g <- runs[[ref_name]]
  
  #Checkpoint: a run denoised by a previous (interrupted) execution is reloaded
  checkpoint <- file.path(path.output, "denoising", g$run[1],
                          paste0("dada_results_", g$error_function[1], ".rds"))
  if(file.exists(checkpoint)){
    message("  [skip] already denoised: ", checkpoint)
    dada_results[[ref_name]] <- readRDS(checkpoint)
    next
  }
  
  res <- tryCatch({
    err_obj <- readRDS(g$learn_error_rds[1])
    getErrors(err_obj)   # fails early if it is not a valid error object
    dd <- dada(g$filtered_fastq, err = err_obj, multithread = n_cpus)
    if(inherits(dd, "dada")) dd <- list(dd)
    names(dd) <- sub("\\.fastq\\.gz$", "", basename(g$filtered_fastq))
    dd
  }, error = function(e){
    message("  [failed] ", ref_name, ": ", conditionMessage(e))
    NULL
  })
  if(is.null(res)){
    failed[[ref_name]] <- data.frame(ref_name = ref_name, run = g$run[1],
                                     error_function = g$error_function[1])
    next
  }
  dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, checkpoint)
  dada_results[[ref_name]] <- res
}
if(length(failed) > 0){
  write.csv(bind_rows(failed), file.path(path.output, "failed_runs.csv"), row.names = FALSE)
}
if(length(dada_results) == 0) stop("No denoising result")
saveRDS(dada_results, file.path(path.output, "dada_results_list.rds"))


##------------------------------------------------------------------------
## 4. Sequence tables (list), combination, chimera removal, taxonomy
##------------------------------------------------------------------------
seq_tables <- list()
meta.list  <- list()
N <- length(dada_results)
n <- 1
for (ref_name in names(dada_results)){
  cat("[", n, "/", N, "] sequence table : ", ref_name, "\n")
  n <- n + 1
  res <- dada_results[[ref_name]]
  
  if (length(res) == 0){
    message("  [skipped] empty result: ", ref_name)
    next
  }
  
  seq.tab <- makeSequenceTable(res)
  
  fastq.names <- rownames(seq.tab)
  if (is.null(fastq.names)){
    fastq.names <- paste0("sample", seq_len(nrow(seq.tab)))
    message("  [no sample names] ", ref_name, " -> ", paste(fastq.names, collapse = ", "))
  }
  
  info <- runs[[ref_name]][1, ]
  sample.names <- paste0(ref_name, "_", fastq.names)
  meta.list[[ref_name]] <- data.frame(
    fastq.files      = fastq.names,
    platform         = info$dataset,
    dataset.prop     = info$nbases,
    subset.replicate = info$seed,
    error.func       = info$error_function,
    sample.out       = sample.names
  )
  rownames(seq.tab) <- sample.names
  seq_tables[[ref_name]] <- seq.tab
}
saveRDS(seq_tables, file.path(path.output, "sequence_tables_list.rds"))

all.meta.data <- bind_rows(meta.list) |> as.data.frame()
rownames(all.meta.data) <- all.meta.data$sample.out
all.meta.data$sample.out <- NULL
all.meta.data$dataset.prop <- factor(all.meta.data$dataset.prop, levels = prop.levels)
write.csv(all.meta.data, file.path(path.output, "all_metadata.csv"))

combined.seq_tables <- mergeSequenceTables(tables = seq_tables)
saveRDS(combined.seq_tables, file.path(path.output, "all_dataset_sequence_table.rds"))

combined.seq_tables.nochim <- removeBimeraDenovo(combined.seq_tables, method = "consensus",
                                                 multithread = n_cpus)
saveRDS(combined.seq_tables.nochim, file.path(path.output, "all_dataset_sequence_table_nochim.rds"))

taxaAssign <- assignTaxonomy(combined.seq_tables.nochim, ref.db, multithread = n_cpus)
saveRDS(taxaAssign, file.path(path.output, "all_dataset_taxonomy_assignment.rds"))

message("Done.")
quit(save = "no", status = 0)