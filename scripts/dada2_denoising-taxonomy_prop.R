#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr")

#Check the number of parameters
args <- commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  print("This script takes the following parameters:")
  print("   (1) Path to the directory containing the run directories (<run>/Filtered/*.fastq.gz and <run>/RDS/learn_error_data_<func>.rds)")
  print("   (2) Path to output directory")
  print("   (3) Path to the taxonomy reference database (DEFAULT SILVA v138.2 trainset)")
  stop("Error: Required arguments are not provided.", call.=FALSE)
}

path.runs   <- args[1]
path.output <- args[2]
ref.db      <- if (length(args) >= 3) args[3] else
  "refSeq/SILVA-v138.2-16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz"

if(!dir.exists(path.runs)) stop("Directory not found: ", path.runs)
if(!file.exists(ref.db))   stop("Reference database not found: ", ref.db)
dir.create(path.output, recursive = TRUE, showWarnings = FALSE)

n_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
if (is.na(n_cpus)) n_cpus <- TRUE

#Constants
prop.levels <- c("1", "5", "10", "25", "50", "100")   # % of reads kept
full.prop   <- "100"   # proportion of the runs without "_<proportion>_<seed>" suffix (full dataset)
## Run directory: <dataset>[_<proportion>_<seed>]
prop.alt    <- paste(setdiff(prop.levels, full.prop), collapse = "|")
prop.pat    <- paste0("^(.+?)_(", prop.alt, ")_([0-9]+)$")


##------------------------------------------------------------------------
## 1. Dataset table : dataset x error function x proportion x seed (x fastq file)
##    Each run directory holds its own cleaned reads (Filtered/) and its learned errors (RDS/)
##------------------------------------------------------------------------
err.paths <- list.files(path.runs, pattern = "^learn_error_data_.*\\.rds$",
                        recursive = TRUE, full.names = TRUE)
if(length(err.paths) == 0) stop("No learn_error_data_*.rds found in ", path.runs)

run.dir    <- dirname(dirname(err.paths))
run        <- basename(run.dir)
has.suffix <- grepl(prop.pat, run, perl = TRUE)
seed       <- ifelse(has.suffix, sub(prop.pat, "\\3", run, perl = TRUE), NA_character_)

runs.table <- data.frame(
  dataset          = ifelse(has.suffix, sub(prop.pat, "\\1", run, perl = TRUE), run),
  error_function   = sub("^learn_error_data_", "", sub("\\.rds$", "", basename(err.paths))),
  proportion       = ifelse(has.suffix, sub(prop.pat, "\\2", run, perl = TRUE), full.prop),
  seed             = seed,
  run              = run,
  run_dir          = normalizePath(run.dir),
  learn_error_rds  = normalizePath(err.paths),
  stringsAsFactors = FALSE
)
runs.table$proportion <- factor(runs.table$proportion, levels = prop.levels)  # numeric order when sorting

## One row per cleaned fastq file: <run>/Filtered/*.fastq.gz
datasets.table <- bind_rows(lapply(seq_len(nrow(runs.table)), function(i){
  files <- list.files(file.path(runs.table$run_dir[i], "Filtered"),
                      pattern = "\\.fastq\\.gz$", full.names = TRUE)
  if(length(files) == 0){
    message("[skipped] no fastq.gz in ", file.path(runs.table$run_dir[i], "Filtered"))
    return(NULL)
  }
  rows <- runs.table[rep(i, length(files)), ]
  rows$filtered_fastq <- normalizePath(files)
  rows
}))
if(nrow(datasets.table) == 0) stop("No Filtered/*.fastq.gz file found in ", path.runs)

datasets.table <- datasets.table |> arrange(dataset, proportion, seed, error_function, filtered_fastq)
rownames(datasets.table) <- NULL
write.csv(datasets.table, file.path(path.output, "datasets_table.csv"), row.names = FALSE)
message(n_distinct(datasets.table$dataset), " datasets, ",
        n_distinct(datasets.table$error_function), " error functions, ",
        n_distinct(datasets.table$proportion), " proportions: ",
        n_distinct(datasets.table[c("run", "error_function")]), " runs to denoise")


##------------------------------------------------------------------------
## 2. Denoising : list of results, one element per dataset x proportion x seed x error function
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
## 3. Sequence tables (list), combination, chimera removal, taxonomy
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
    dataset.prop     = as.character(info$proportion),
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