#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr") 
library(wavethresh, verbose = FALSE, quietly = TRUE); packageVersion("wavethresh")

#Constants
## Number of bases used for learning errors. "1e+08" = default (full dataset), 
## i.e. runs whose directory name has no "_1e+0X[_<seed>]" suffix.
prop.levels <- c("1e+04", "1e+05", "1e+06", "1e+07", "1e+08")
full.prop <- "1e+08"
## Run directory: <dataset>[_<nbases>[_<seed>]]  (several seeds per nbases)
nbases.pat <- "^(.+?)_(1e\\+0[0-9])(?:_([0-9]+))?$"
## Runs of the proportion experiment (<dataset>_<prop>_<seed>), not part of this analysis
prop.pat <- "^(.+?)_(1|5|10|25|50)_([0-9]+)$"
## Resource metrics written by track_resources() in dada2_analysis.R
resource.cols <- c("wall_sec", "cpu_sec", "avg_cores_used", "threads_available",
                   "cpus_allocated", "peak_ram_gb", "end_ram_gb")

#Define functions

## Parse <dataset>[_<nbases>[_<seed>]]/RDS/<prefix>_<errors.function>.<ext>
## Dataset names contain underscores (e.g. LIB_16S_KINNEX_SEGMENTED_REVIO_SPRQ_NX), 
## so we cannot split on "_": the nbases (and optional seed) suffix is matched 
## from the end of the run directory name.
parseRunPath <- function(file_path){
  run  <- basename(dirname(dirname(file_path)))
  file <- basename(file_path)
  has.suffix <- grepl(nbases.pat, run, perl = TRUE)
  
  used.seed <- ifelse(has.suffix, sub(nbases.pat, "\\3", run, perl = TRUE), NA_character_)
  used.seed[!is.na(used.seed) & used.seed == ""] <- NA_character_
  
  errors.function <- sub("^(dada_results|learn_error_data|execution_time)_", "", file)
  errors.function <- sub("\\.(rds|csv)$", "", errors.function)
  
  data.frame(
    ref_name        = paste0(run, "_", file),
    run             = run,
    dataset         = ifelse(has.suffix, sub(nbases.pat, "\\1", run, perl = TRUE), run),
    dataset.prop    = ifelse(has.suffix, sub(nbases.pat, "\\2", run, perl = TRUE), full.prop),
    used.seed       = used.seed,
    errors.function = errors.function,
    stringsAsFactors = FALSE
  )
}

readMultiRDS <- function(list_of_paths){
  res <- list()
  for(file_path in list_of_paths){
    ref_name <- paste0(basename(dirname(dirname(file_path))), "_", basename(file_path))
    res[[ref_name]] <- readRDS(file_path)
  }
  return(res)
}

## completed: ref_names of the dada_results RDS that exist. A run without a 
## dada_results file (e.g. loessErrfun_mod1-4, makeBinnedQualErrfun on Sequel) 
## only logged the failed learnError step, so it is flagged as "failed".
readMultiCSV <- function(list_of_paths, completed = character()){
  res <- list()
  for(file_path in list_of_paths){
    info <- parseRunPath(file_path)
    df <- read.csv(file_path)
    if(all(df$process_name == "learnError")) next
    df$X <- NULL
    df$platform          <- info$dataset
    df$dataset.prop      <- info$dataset.prop
    df$used.seed         <- info$used.seed
    df$errors.function   <- info$errors.function
    df$status            <- ifelse(paste0(info$run, "_dada_results_", info$errors.function, ".rds") %in% completed,
                                   "completed", "failed")
    res[[info$ref_name]] <- df
  }
  res <- bind_rows(res)
  
  #Older runs only have process_name/start_time/end_time
  for(col in setdiff(resource.cols, names(res))) res[[col]] <- NA_real_
  res[resource.cols] <- lapply(res[resource.cols], as.numeric)
  
  res$dataset.prop <- factor(res$dataset.prop, levels = prop.levels)
  res$end_time <- as.POSIXct(res$end_time)
  res$start_time <- as.POSIXct(res$start_time)
  #Wall time (s): use the logged value, fall back to end - start
  res$duration <- ifelse(!is.na(res$wall_sec), res$wall_sec,
                         as.numeric(difftime(res$end_time, res$start_time, units = "secs")))
  #Cores really allocated (SLURM) if known, otherwise the threads visible to RcppParallel
  res$cores_ref      <- coalesce(res$cpus_allocated, res$threads_available)
  res$cpu_efficiency <- res$avg_cores_used / res$cores_ref
  res$status <- factor(res$status, levels = c("completed", "failed"))
  return(
    res
  )
}

## Distance between the error matrix learned on the full dataset and the one 
## learned on a subset of bases (one value per seed), for every error function 
## available in the dataset.
## Returns a wide data.frame (one column "<errors.function>.rds" per function), 
## one row per (dataset.prop, used.seed).
safeGetErrors <- function(dada_result, ref_name){
  tryCatch(getErrors(dada_result),
           error = function(e){
             message("  [skipped] ", ref_name, ": ", conditionMessage(e))
             NULL
           })
}

computeDistances <- function(dada2_results_data, dada2_meta, dataset_name){
  meta <- dada2_meta[dada2_meta$dataset == dataset_name, ]
  distances <- list()
  for(err.func in unique(meta$errors.function)){
    sub.meta <- meta[meta$errors.function == err.func, ]
    full.ref <- sub.meta$ref_name[sub.meta$dataset.prop == full.prop]
    if(length(full.ref) != 1){
      message("[", dataset_name, "] no full-dataset result for ", err.func, ": skipped")
      next
    }
    full <- safeGetErrors(dada2_results_data[[full.ref]], full.ref)
    if(is.null(full)){
      message("[", dataset_name, "] full-dataset error matrix is NULL for ", err.func, ": skipped")
      next
    }
    for(i in seq_len(nrow(sub.meta))){
      other <- safeGetErrors(dada2_results_data[[sub.meta$ref_name[i]]], sub.meta$ref_name[i])
      if(is.null(other)) next
      if(!identical(dim(full), dim(other))){
        message("[", dataset_name, "] dimension mismatch for ", sub.meta$ref_name[i], ": skipped")
        next
      }
      distances[[length(distances) + 1]] <- data.frame(
        dataset.prop    = sub.meta$dataset.prop[i],
        used.seed       = sub.meta$used.seed[i],
        errors.function = paste0(err.func, ".rds"),
        distances       = l2norm(full, other),   # L2 norm
        stringsAsFactors = FALSE
      )
    }
  }
  if(length(distances) == 0) return(NULL)
  distances <- bind_rows(distances) |>
    pivot_wider(names_from = errors.function, values_from = distances) |>
    as.data.frame()
  #dataset.prop stays a character ("1e+04", ...) to match prop.levels exactly
  distances$used.seed <- as.numeric(distances$used.seed)
  return(distances)
}


transToLong.prop <- function(distances){
  distances.long <- distances |>
    pivot_longer(cols = ends_with(".rds"),
                 values_to = "distances",
                 names_to = "errors.function",
                 values_drop_na = TRUE) |>
    mutate(
      errors.function = sub(".rds","",errors.function, fixed = TRUE),
      dataset.prop = factor(dataset.prop, levels = prop.levels)
    )
  return(distances.long)
}

plotDistancePoint.prop <- function(distances.long){
  distances.long |>
    ggplot()+
    geom_point(aes(x=dataset.prop, y=distances), alpha=0.4)+
    facet_grid(cols=vars(errors.function))+
    labs(x="Number of bases for learning errors",
         y="Distance to the full dataset") + theme_bw()
}

plotDistanceBoxplot.prop <- function(distances.long){
  distances.long |>
    ggplot()+
    geom_boxplot(aes(x=dataset.prop, y=distances), alpha=0.4)+
    facet_grid(cols=vars(errors.function))+
    labs(x="Number of bases for learning errors",
         y="Distance to the full dataset") + theme_bw()
}

## Generic resource plot: one boxplot per nbases (spread over the seeds), 
## individual runs as points (shape = HPC node), median as line
plotResourcePlot.prop <- function(execution.time.long, metric, ylab){
  dodge <- position_dodge(width = 0.75)
  execution.time.long |>
    ggplot(aes(x = dataset.prop, y = .data[[metric]], colour = process_name))+
    geom_boxplot(alpha = 0.6, outliers = FALSE, position = dodge)+
    geom_point(aes( group = process_name), position = dodge, alpha = 0.5)+
    stat_summary(aes(group = process_name), fun = median, geom = "line", position = dodge)+
    facet_grid(cols=vars(errors.function), rows = vars(label),
               labeller = labeller(errors.function = label_wrap_gen(width = 100),
                                   label = label_wrap_gen(width = 20)))+
    labs(x="Number of bases for learning errors",
         y=ylab,
         colour = "Process name") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 0))
}

plotExecutionTimePlot.prop <- function(execution.time.long){
  plotResourcePlot.prop(execution.time.long, "duration", "Execution Time (s)")
}

## Size the figure from the number of facets (many datasets now)
saveFacetPlot <- function(plot, filename, data, by.platform = TRUE, ...){
  n.col <- n_distinct(data$errors.function)
  n.row <- if(by.platform) n_distinct(data$platform) else 3
  ggsave(file.path(results.path, filename), plot = plot,
         width = 2.5 * n.col + 2, height = 2 * n.row + 1, limitsize = FALSE, ...)
}

## Is the distance to the full dataset significantly > 0 at each number of bases?
## (t-test if Shapiro-Wilk does not reject normality, Wilcoxon otherwise), 
## using the seeds as replicates. Needs >= 3 seeds per nbases (Shapiro-Wilk); 
## groups with fewer are dropped. Returns NULL if there is nothing to test.
runDistanceTests <- function(distances){
  results <- list()
  for(err.func in grep("\\.rds$", names(distances), value = TRUE)){
    test_by_group <- distances |>
      #The full dataset is at distance 0 from itself: nothing to test
      filter(dataset.prop != full.prop, !is.na(.data[[err.func]])) |>
      mutate(dataset.prop = factor(dataset.prop, levels = setdiff(prop.levels, full.prop))) |>
      group_by(dataset.prop) |>
      summarise(
        n         = n(),
        mean      = mean(.data[[err.func]], na.rm = TRUE),
        sd        = sd(.data[[err.func]], na.rm = TRUE),
        median    = median(.data[[err.func]], na.rm = TRUE),
        shapiro.p = tryCatch(shapiro.test(.data[[err.func]])$p.value,
                             error = function(e) NA_real_),
        p.value   = if_else(shapiro.p > 0.05, tryCatch(t.test(.data[[err.func]], mu = 0, alternative = "greater")$p.value,
                                                       error = function(e) NA_real_),
                            tryCatch(wilcox.test(.data[[err.func]], mu = 0, alternative = "greater")$p.value,
                                     error = function(e) NA_real_)),
        .groups = "drop"
      )|>
      mutate(err.func = err.func)
    
    results[[err.func]] <- test_by_group
  }
  if(length(results) == 0) return(NULL)
  all_results <- bind_rows(results) |>
    mutate(p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
           p.adj.benjamin.hochberg = p.adjust(p.value, method = "BH"))
  all_results <- all_results |> drop_na(shapiro.p)
  return(all_results)
}

#Input data
path.alldataset <- file.path(list.files("results", full.names = T), "RDS")
paths_learn_errors_data <- list.files(path.alldataset, pattern = "learn_error", full.names = T) 
paths_dada2_results_data <- list.files(path.alldataset, pattern = "dada_results", full.names = T)
paths_exec_times <- list.files(path.alldataset, pattern = "execution_time", full.names = T) 
ref.db <- file.path("refSeq/SILVA-v138.2-16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz")


datasets_details <- read.csv("summary/Dataset_Summary-Before_Normalization.csv")
datasets_details <- datasets_details |>
  rowwise() |>
  mutate(
    label = {
      optional <- c(Library.Type, Segmentation, Chemistry)
      optional <- optional[!is.na(optional) & optional != ""]
      paste0(Location, " (", paste(c(Sequencer, optional), collapse = "/"), ")")
    }
  ) |>
  ungroup()

dataset.labels <- setNames(datasets_details$label, datasets_details$Dataset.Name)

getLabel <- function(x){
  out <- unname(dataset.labels[x])
  ifelse(is.na(out), x, out)   # fall back to the raw name if no match
}

## Learning-error logs: keep the nbases runs only, and parse the run name 
## (<dataset>[_<nbases>[_<seed>]]) so that every seed is a separate row.
learn_error_logs <- read.csv("summary/new_results_summary_learning.csv")
learn_error_logs <- learn_error_logs[!grepl(prop.pat, learn_error_logs$Dataset, perl = TRUE),]
learn_error_logs$nbases[is.na(learn_error_logs$nbases)] <- 1e+08
learn_error_logs <- learn_error_logs |>
  mutate(
    has.suffix = grepl(nbases.pat, Dataset, perl = TRUE),
    platform   = ifelse(has.suffix, sub(nbases.pat, "\\1", Dataset, perl = TRUE), Dataset),
    used.seed  = ifelse(has.suffix, sub(nbases.pat, "\\3", Dataset, perl = TRUE), NA_character_),
    nbases     = factor(sprintf("%.0e", nbases), levels = prop.levels),
    label      = getLabel(platform)
  )
learn_error_logs$used.seed[learn_error_logs$used.seed %in% ""] <- NA_character_

results.path <- "summary/nbases"
dir.create(results.path, showWarnings = F)

#Load data
learn_errors_data <- readMultiRDS(paths_learn_errors_data)
dada2_results_data <- readMultiRDS(paths_dada2_results_data)
exec_times_data <- readMultiCSV(paths_exec_times, completed = names(dada2_results_data))

#Run metadata (dataset / nbases / seed / error function), indexed by ref_name
dada2_meta <- bind_rows(lapply(paths_dada2_results_data, parseRunPath)) |> as.data.frame()
rownames(dada2_meta) <- dada2_meta$ref_name
datasets <- unique(dada2_meta$dataset)
message(length(datasets), " datasets found:\n  ", paste(datasets, collapse = "\n  "))

#Number of seeds (replicates) available for each dataset / nbases / error function
seed_counts <- dada2_meta |>
  group_by(dataset, dataset.prop, errors.function) |>
  summarise(n_seeds = n_distinct(used.seed, na.rm = TRUE), .groups = "drop")
write.csv(seed_counts, file.path(results.path, "nbases_seed_counts.csv"), row.names = F)

#Failed runs (no dada_results file: only the failed learnError step was logged)
failed_runs <- exec_times_data |>
  filter(status == "failed") |>
  distinct(platform, dataset.prop, used.seed, errors.function)
write.csv(failed_runs, file.path(results.path, "nbases_failed_runs.csv"), row.names = F)
failed_runs |> count(errors.function, name = "n_failed_runs") |> print()

#Execution time and resource usage (completed runs only)
exec_times_ok <- exec_times_data |> filter(status == "completed")
exec_times_ok$label <- getLabel(exec_times_ok$platform)
#One run = dataset + nbases + seed + error function
exec_times_ok <- exec_times_ok |> mutate(
  id = paste(platform, dataset.prop, used.seed, errors.function, sep = "_")
)

tmp <- learn_error_logs |>
  mutate(id = paste(platform, nbases, used.seed, Loss.Function, sep = "_")) |>
  dplyr::select(id, HPC.node)
if(anyDuplicated(tmp$id)){
  warning("Several learning-error log rows share the same dataset/nbases/seed/function: ",
          "only the first HPC node is kept for each.")
  tmp <- distinct(tmp, id, .keep_all = TRUE)
}

exec_times_ok <- left_join(exec_times_ok, tmp, by=c("id"="id")) 
exec_times_ok$HPC.node[is.na(exec_times_ok$HPC.node)] <- "unknown"

## Execution time
p <- plotExecutionTimePlot.prop(exec_times_ok)
ggsave(file.path(results.path, "nbases_Execution_Time_boxplot.svg"))
#saveFacetPlot(p, "nbases_Execution_Time_boxplot.pdf", exec_times_ok)

## Total CPU time
p <- plotResourcePlot.prop(exec_times_ok, "cpu_sec", "CPU Time (s)")
ggsave(file.path(results.path, "nbases_CPU_Time_plot.svg"))
#saveFacetPlot(p, "nbases_CPU_Time_plot.pdf", exec_times_ok)

## Average number of busy cores (cpu_sec / wall_sec)
p <- plotResourcePlot.prop(exec_times_ok, "avg_cores_used", "Average cores used")
ggsave(file.path(results.path, "nbases_Cores_Used_plot.svg"))
#saveFacetPlot(p, "nbases_Cores_Used_plot.pdf", exec_times_ok)

## CPU efficiency: average busy cores / cores allocated
p <- plotResourcePlot.prop(exec_times_ok, "cpu_efficiency", "CPU efficiency (avg cores used / cores available)")
ggsave(file.path(results.path, "nbases_CPU_Efficiency_plot.svg"))
#saveFacetPlot(p, "nbases_CPU_Efficiency_plot.pdf", exec_times_ok)

## Peak RAM
## NB: peak is per step only if /proc/self/clear_refs was writable; otherwise 
## VmHWM is cumulative and Denoising peak >= learnError peak.
p <- plotResourcePlot.prop(exec_times_ok, "peak_ram_gb", "Peak RAM (GB)")
ggsave(file.path(results.path, "nbases_Peak_RAM_plot.svg"))
#saveFacetPlot(p, "nbases_Peak_RAM_plot.pdf", exec_times_ok)

## Summary table (mean, sd and median over replicates/seeds)
resource_summary <- exec_times_ok |>
  group_by(platform, errors.function, dataset.prop, process_name) |>
  summarise(
    n_runs = n(),
    across(c(duration, cpu_sec, avg_cores_used, cpu_efficiency, peak_ram_gb, end_ram_gb),
           list(mean   = ~ mean(.x, na.rm = TRUE),
                sd     = ~ sd(.x, na.rm = TRUE),
                median = ~ median(.x, na.rm = TRUE)),
           .names = "{.fn}_{.col}"),
    .groups = "drop"
  )
write.csv(resource_summary, file.path(results.path, "nbases_resource_summary.csv"), row.names = F)

#Distances to the full dataset, for every dataset (one value per seed)
distances_by_dataset <- list()
distances_long_by_dataset <- list()
for(ds in datasets){
  ds.distances <- computeDistances(dada2_results_data, dada2_meta, ds)
  if(is.null(ds.distances)) next
  ds.distances.long <- transToLong.prop(ds.distances)
  distances_by_dataset[[ds]] <- ds.distances
  distances_long_by_dataset[[ds]] <- ds.distances.long |> mutate(platform = ds)
  
  #p <- plotDistancePoint.prop(ds.distances.long)
  #ggsave(file.path(results.path, paste0("nbases_", ds, "_dataset_distance_plot.pdf")), plot = p)
  #p <- plotDistanceBoxplot.prop(ds.distances.long)
  #ggsave(file.path(results.path, paste0("nbases_", ds, "_dataset_distance_boxplot.pdf")), plot = p)
}

#Combination of datasets
combined_distances.long <- bind_rows(distances_long_by_dataset)
write.csv(combined_distances.long, file.path(results.path, "nbases_all_distances.csv"), row.names = F)

combined_distances.long$label <- getLabel(combined_distances.long$platform)

p <- combined_distances.long |>
  ggplot(aes(x = dataset.prop, y = distances)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(colour = dataset.prop), width = 0.15, alpha = 0.6) +
  stat_summary(aes(group = 1), fun = median, geom = "line", colour = "grey30") +
  facet_grid(cols = vars(errors.function), rows = vars(label),
             labeller = labeller(errors.function = label_wrap_gen(width = 100),
                                 label = label_wrap_gen(width = 20))) +
  labs(x = "Number of bases for learning errors",
       y = "Distance to the full dataset") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = "none")
ggsave(file.path(results.path, "Combined_dataset_distance_boxplot.svg"))
#saveFacetPlot(p, "Combined_dataset_distance_boxplot.pdf", combined_distances.long)

#Statistics tests
## One test table per dataset (written as nbases_<dataset>_summary_test.txt, lower case)
test_results <- list()
for(ds in names(distances_by_dataset)){
  ds.results <- runDistanceTests(distances_by_dataset[[ds]])
  if(is.null(ds.results) || nrow(ds.results) == 0){
    message("[", ds, "] not enough replicates (seeds) for statistical tests: skipped")
    next
  }
  write.table(ds.results, file.path(results.path, paste0("nbases_", tolower(ds), "_summary_test.txt")), row.names = F)
  ds.results |>
    dplyr::filter(p.adj.bonferroni > 0.01 | p.adj.benjamin.hochberg > 0.01) |>
    print()
  test_results[[ds]] <- ds.results
}

##Combination and visualization 
if(length(test_results) > 0){
  all_platforms <- bind_rows(test_results, .id = "platform") |>
    mutate(
      dataset.prop = factor(dataset.prop, levels = setdiff(prop.levels, full.prop)),
      signif = case_when(
        p.adj.bonferroni < 0.001 ~ "< 0.001",
        p.adj.bonferroni < 0.01  ~ "< 0.01",
        p.adj.bonferroni < 0.05  ~ "< 0.05",
        TRUE ~ "ns"
      )
    )
  all_platforms$label <- getLabel(all_platforms$platform)
  
  p <- ggplot(all_platforms, aes(x = dataset.prop, y = mean, group = err.func, colour = err.func)) +
    geom_line(position = position_dodge(width = 0.3)) +
    geom_point(aes(shape = signif), size = 4, position = position_dodge(width = 0.3)) +
    facet_wrap(~label, scales = "free_y", ncol = 3,
               labeller = labeller(label = label_wrap_gen(width = 30)))  +
    labs(x = "Number of bases for learning errors", y = "Average distance to the full dataset",
         color = "Error Function", shape = "Adjusted p.value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(results.path, "nbases_summary_statistics_test.pdf"), plot = p,
         width = 14, height = 4 * ceiling(n_distinct(all_platforms$platform) / 3), limitsize = FALSE)
}

#Plot learning Error log data (one point per seed)
ggplot(learn_error_logs,
       aes(x = nbases, y = Total.bases.used, colour = label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  facet_wrap(~Loss.Function, ncol = 4) +
  labs(x = "Number of bases for learning errors", y = "Total bases used", colour = "Dataset") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(results.path, "Total_bases_used_during_learingError.pdf"))

ggplot(learn_error_logs,
       aes(x = nbases, y = Total.reads.used, colour = label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  facet_wrap(~Loss.Function, ncol = 4) +
  labs(x = "Number of bases for learning errors", y = "Total number of reads used", colour = "Dataset") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(results.path, "Total_reads_used_during_learingError.pdf"))

ggplot(learn_error_logs, aes(x = Total.bases.used, y = label)) +
  geom_point(aes(size  = Total.reads.used, colour = nbases), alpha = 0.5) +
  scale_size_continuous(range = c(0.5, 5)) +
  facet_wrap(~Loss.Function, ncol = 4) +
  labs(size = "Total reads used", y = NULL,
       x="Total bases used", 
       colour = "Number of bases for learning errors") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(results.path, "Summary_learingError_reads-bases.pdf"))

#Sequence table
seq_tables <- list()
meta.list  <- list()
N <- length(names(dada2_results_data))
n <- 1
for (ref_name in names(dada2_results_data)){
  cat("[", n, "/", N, "] : ", ref_name, "\n")
  n <- n + 1
  res <- dada2_results_data[[ref_name]]
  
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
  
  info <- dada2_meta[ref_name, ]
  sample.names <- paste0(sub(".rds", "", ref_name, fixed = TRUE), "_", fastq.names)
  meta.list[[ref_name]] <- data.frame(
    fastq.files      = fastq.names,
    platform         = info$dataset,
    dataset.prop     = info$dataset.prop,
    subset.replicate = info$used.seed,
    error.func       = info$errors.function,
    sample.out       = sample.names
  )
  rownames(seq.tab) <- sample.names
  seq_tables[[ref_name]] <- seq.tab
}
all.meta.data <- bind_rows(meta.list) |> as.data.frame()
rownames(all.meta.data) <- all.meta.data$sample.out
all.meta.data$sample.out <- NULL
all.meta.data$dataset.prop <- factor(all.meta.data$dataset.prop, levels = prop.levels)
write.csv(all.meta.data,file.path(results.path, "nbases_all_metadata.csv"))
combined.seq_tables <- mergeSequenceTables(tables=seq_tables)
saveRDS(combined.seq_tables,file.path(results.path, "nbases_all_dataset_sequence_table.rds"))

#Assign taxonomy
combined.seq_tables.nochim <- removeBimeraDenovo(combined.seq_tables, method="consensus", multithread=TRUE)
saveRDS(combined.seq_tables.nochim,file.path(results.path, "nbases_all_dataset_sequence_table_nochim.rds"))
taxaAssign <- assignTaxonomy(combined.seq_tables.nochim, ref.db,
                             multithread = T)
saveRDS(taxaAssign,file.path(results.path, "nbases_all_dataset_taxonomy_assignment.rds"))

message("Done.")
quit(save = "no", status = 0)