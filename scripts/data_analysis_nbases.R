
#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr") 

#Constants
## Number of bases used for learning errors. "1e+08" = default (full dataset), 
## i.e. runs whose directory name has no "_1e+0X" suffix.
prop.levels <- c("1e+04", "1e+05", "1e+06", "1e+07", "1e+08")
## Resource metrics written by track_resources() in dada2_analysis.R
resource.cols <- c("wall_sec", "cpu_sec", "avg_cores_used", "threads_available",
                   "cpus_allocated", "peak_ram_gb", "end_ram_gb")

#Define functions

## Parse <dataset>[_<nbases>][_<seed>]/RDS/<prefix>_<errors.function>.<ext>
## Dataset names now contain underscores (e.g. LIB_16S_KINNEX_SEGMENTED_REVIO_SPRQ_NX), 
## so we cannot split on "_" anymore: the nbases (and optional seed) suffix is matched 
## from the end of the run directory name.
parseRunPath <- function(file_path){
  run  <- basename(dirname(dirname(file_path)))
  file <- basename(file_path)
  pat  <- "^(.+?)_(1e\\+0[0-9])(?:_([0-9]+))?$"
  has.suffix <- grepl(pat, run, perl = TRUE)
  
  used.seed <- ifelse(has.suffix, sub(pat, "\\3", run, perl = TRUE), NA_character_)
  used.seed[!is.na(used.seed) & used.seed == ""] <- NA_character_
  
  errors.function <- sub("^(dada_results|learn_error_data|execution_time)_", "", file)
  errors.function <- sub("\\.(rds|csv)$", "", errors.function)
  
  data.frame(
    ref_name        = paste0(run, "_", file),
    run             = run,
    dataset         = ifelse(has.suffix, sub(pat, "\\1", run, perl = TRUE), run),
    dataset.prop    = ifelse(has.suffix, sub(pat, "\\2", run, perl = TRUE), "1e+08"),
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
    df$platform        <- info$dataset
    df$dataset.prop    <- info$dataset.prop
    df$used.seed       <- info$used.seed
    df$errors.function <- info$errors.function
    df$status <- ifelse(paste0(info$run, "_dada_results_", info$errors.function, ".rds") %in% completed,
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
## learned on a subset of bases, for every error function available in the dataset.
## Returns a wide data.frame (one column "<errors.function>.rds" per function).
computeDistances <- function(dada2_results_data, dada2_meta, dataset_name){
  meta <- dada2_meta[dada2_meta$dataset == dataset_name, ]
  distances <- list()
  for(err.func in unique(meta$errors.function)){
    sub.meta <- meta[meta$errors.function == err.func, ]
    full.ref <- sub.meta$ref_name[sub.meta$dataset.prop == "1e+08"]
    if(length(full.ref) != 1){
      message("[", dataset_name, "] no full-dataset result for ", err.func, ": skipped")
      next
    }
    full <- getErrors(dada2_results_data[[full.ref]])
    for(i in seq_len(nrow(sub.meta))){
      other <- getErrors(dada2_results_data[[sub.meta$ref_name[i]]])
      if(!identical(dim(full), dim(other))){
        message("[", dataset_name, "] dimension mismatch for ", sub.meta$ref_name[i], ": skipped")
        next
      }
      distances[[length(distances) + 1]] <- data.frame(
        dataset.prop    = sub.meta$dataset.prop[i],
        used.seed       = sub.meta$used.seed[i],
        errors.function = paste0(err.func, ".rds"),
        distances       = sqrt(sum((full - other)^2)),   # L2 norm
        stringsAsFactors = FALSE
      )
    }
  }
  if(length(distances) == 0) return(NULL)
  distances <- bind_rows(distances) |>
    pivot_wider(names_from = errors.function, values_from = distances) |>
    as.data.frame()
  distances$dataset.prop <- as.numeric(distances$dataset.prop)
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
      dataset.prop = factor(format(dataset.prop, scientific = TRUE), 
                            levels = prop.levels)
    )
  return(distances.long)
}

plotDistancePoint.prop <- function(distances.long){
  distances.long |>
    ggplot(aes(x=dataset.prop, y=distances))+
    geom_jitter(alpha=0.6, size = 3)+
    facet_grid(cols=vars(errors.function))+
    labs(x="Number bases to consider for learning error",
         y="Distance to the default number") + theme_bw()
}

plotDistanceBoxplot.prop <- function(distances.long){
  distances.long |>
    ggplot(aes(x=dataset.prop, y=distances))+
    geom_boxplot(alpha=0.4, outliers = F)+
    geom_jitter(alpha=0.6, size = 3)+
    facet_grid(cols=vars(errors.function))+
    labs(x="Number bases to consider for learning error",
         y="Distance to the default number") + theme_bw()
}

## Generic resource plot: one point per run/step, metric on y
plotResourcePlot.prop <- function(execution.time.long, metric, ylab){
  dodge <- position_dodge(width = 0.3)
  n.platform <- n_distinct(execution.time.long$platform)
  execution.time.long |>
    ggplot(aes(x=dataset.prop, y=.data[[metric]], colour = process_name,
               group = interaction(platform, process_name)))+
    geom_line(position = dodge, alpha=0.4)+
    geom_point(aes(shape=platform), position = dodge, alpha=0.6)+
    scale_shape_manual(values = rep_len(c(16, 17, 15, 3, 7, 8, 4, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14), n.platform))+
    facet_grid(cols=vars(errors.function), scales = "free",
               labeller = labeller(errors.function = label_wrap_gen(width = 100)))+
    labs(x="Number of bases for learning errors",
         y=ylab,
         colour = "Process Name",
         shape = "Dataset") + theme_bw()
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

#Input data
path.alldataset <- file.path(list.files("results_nbases", full.names = T), "RDS")
paths_learn_errors_data <- list.files(path.alldataset, pattern = "learn_error", full.names = T) 
paths_dada2_results_data <- list.files(path.alldataset, pattern = "dada_results", full.names = T)
paths_exec_times <- list.files(path.alldataset, pattern = "execution_time", full.names = T) 
ref.db <- file.path("refSeq/SILVA-v138.2-16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz")

results.path <- "summary"
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

#Failed runs (no dada_results file: only the failed learnError step was logged)
failed_runs <- exec_times_data |>
  filter(status == "failed") |>
  distinct(platform, dataset.prop, used.seed, errors.function)
write.csv(failed_runs, file.path(results.path, "nbases_failed_runs.csv"), row.names = F)
failed_runs |> count(errors.function, name = "n_failed_runs") |> print()

#Execution time and resource usage (completed runs only)
exec_times_ok <- exec_times_data |> filter(status == "completed")

## Execution time
p <- plotExecutionTimePlot.prop(exec_times_ok)
saveFacetPlot(p, "nbases_Execution_Time_boxplot.pdf", exec_times_ok)

## Total CPU time
p <- plotResourcePlot.prop(exec_times_ok, "cpu_sec", "CPU Time (s)")
saveFacetPlot(p, "nbases_CPU_Time_plot.pdf", exec_times_ok)

## Average number of busy cores (cpu_sec / wall_sec)
p <- plotResourcePlot.prop(exec_times_ok, "avg_cores_used", "Average cores used")
saveFacetPlot(p, "nbases_Cores_Used_plot.pdf", exec_times_ok)

## CPU efficiency: average busy cores / cores allocated
p <- plotResourcePlot.prop(exec_times_ok, "cpu_efficiency", "CPU efficiency (avg cores used / cores available)")
saveFacetPlot(p, "nbases_CPU_Efficiency_plot.pdf", exec_times_ok)

## Peak RAM
## NB: peak is per step only if /proc/self/clear_refs was writable; otherwise 
## VmHWM is cumulative and Denoising peak >= learnError peak.
p <- plotResourcePlot.prop(exec_times_ok, "peak_ram_gb", "Peak RAM (GB)")
saveFacetPlot(p, "nbases_Peak_RAM_plot.pdf", exec_times_ok)

## Summary table (mean over replicates/seeds)
resource_summary <- exec_times_ok |>
  group_by(platform, errors.function, dataset.prop, process_name) |>
  summarise(
    n_runs = n(),
    across(c(duration, cpu_sec, avg_cores_used, cpu_efficiency, peak_ram_gb, end_ram_gb),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    .groups = "drop"
  )
write.csv(resource_summary, file.path(results.path, "nbases_resource_summary.csv"), row.names = F)

#Distances to the full dataset, for every dataset
distances_by_dataset <- list()
distances_long_by_dataset <- list()
for(ds in datasets){
  ds.distances <- computeDistances(dada2_results_data, dada2_meta, ds)
  if(is.null(ds.distances)) next
  ds.distances.long <- transToLong.prop(ds.distances)
  distances_by_dataset[[ds]] <- ds.distances
  distances_long_by_dataset[[ds]] <- ds.distances.long |> mutate(platform = ds)
  
  p <- plotDistancePoint.prop(ds.distances.long)
  ggsave(file.path(results.path, paste0("nbases_", ds, "_dataset_distance_plot.pdf")), plot = p)
  p <- plotDistanceBoxplot.prop(ds.distances.long)
  ggsave(file.path(results.path, paste0("nbases_", ds, "_dataset_distance_boxplot.pdf")), plot = p)
}

#Combination of datasets
combined_distances.long <- bind_rows(distances_long_by_dataset)
write.csv(combined_distances.long, file.path(results.path, "nbases_all_distances.csv"), row.names = F)

p <- combined_distances.long |>
  ggplot()+
  #geom_boxplot(aes(x=dataset.prop, y=distances), outliers = F)+
  geom_jitter(aes(x=dataset.prop, y=distances, colour = dataset.prop), alpha=0.6)+
  facet_grid(cols=vars(errors.function), rows = vars(platform), scales = "free",
             labeller = labeller(errors.function = label_wrap_gen(width = 100),
                                 platform = label_wrap_gen(width = 20)))+
  labs(x="Number of bases for learning errors",
       y="Distance to the full dataset") + 
  labs(colour="Number of\nbases") + theme_bw()
saveFacetPlot(p, "nbases_Combined_dataset_distance_boxplot.png", combined_distances.long)

dodge <- position_dodge(width = 0.3)
p <- combined_distances.long |>
  ggplot(aes(x=dataset.prop, y=distances, colour = platform, shape = platform, group = platform))+
  geom_line(position = dodge, alpha=0.4)+
  geom_point(position = dodge, alpha=0.6)+
  scale_shape_manual(values = rep_len(c(16, 17, 15, 3, 7, 8, 4, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14),
                                      n_distinct(combined_distances.long$platform)))+
  facet_grid(cols=vars(errors.function), scales = "free",
             labeller = labeller(errors.function = label_wrap_gen(width = 100)))+
  labs(x="Number of bases for learning errors",
       y="Distance to the full dataset",
       colour="Dataset", shape="Dataset")+theme_bw()
saveFacetPlot(p, "nbases_Combined_dataset_distance_plot.pdf", combined_distances.long)


p <- combined_distances.long |>
  dplyr::filter(platform!="Sequel_UniBe") |> 
  ggplot(aes(x=dataset.prop, y=distances, colour = platform, shape = platform, group = platform))+
  geom_line(position = dodge, alpha=0.4)+
  geom_point(position = dodge, alpha=0.6)+
  scale_shape_manual(values = rep_len(c(16, 17, 15, 3, 7, 8, 4, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14),
                                      n_distinct(combined_distances.long$platform)))+
  facet_grid(cols=vars(errors.function), scales = "free",
             labeller = labeller(errors.function = label_wrap_gen(width = 100)))+
  labs(x="Number of bases for learning errors",
       y="Distance to the full dataset",
       colour="Dataset", shape="Dataset")+theme_bw()
saveFacetPlot(p, "nbases_Combined_dataset_distance_plot_without_Sequel_UniBe.pdf", combined_distances.long)

#Statistics tests

# ##Revio_UniBe
# results <- list()
# for(err.func in c("loessErrfun.rds", "PacBioErrfun.rds", 
#                   "makeBinnedQualErrfun.rds", "loessErrfun_mod0.rds")){
  
#   test_by_group <- revio_unibe_distances |>
#     mutate(dataset.prop = factor(format(dataset.prop, scientific = TRUE), levels = c("1e+04", "1e+05", "1e+06", "1e+07"))) |>
#     group_by(dataset.prop) |>
#     summarise(
#       n         = n(),
#       mean      = mean(.data[[err.func]], na.rm = TRUE),
#       shapiro.p = tryCatch(shapiro.test(.data[[err.func]])$p.value,
#                            error = function(e) NA_real_),
#       p.value   = if_else(shapiro.p > 0.05, tryCatch(t.test(.data[[err.func]], mu = 0, alternative = "greater")$p.value,
#                                                      error = function(e) NA_real_),
#                           tryCatch(wilcox.test(.data[[err.func]], mu = 0, alternative = "greater")$p.value,
#                                    error = function(e) NA_real_)),
#       .groups = "drop"
#     )|>
#     mutate(err.func = err.func)
  
#   results[[err.func]] <- test_by_group
# }

# all_results <- bind_rows(results) |>
#   mutate(p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
#          p.adj.benjamin.hochberg = p.adjust(p.value, method = "BH"))
# all_results <- all_results |> drop_na(shapiro.p)
# write.table(all_results, file.path(results.path, "nbases_revio_unibe_summary_test.txt"), row.names = F)

# all_results |>
#   dplyr::filter(p.adj.bonferroni > 0.01 | p.adj.benjamin.hochberg >0.01) |>
#   print()

# ##Sequel_UniBe
# results <- list()
# for(err.func in c("loessErrfun.rds", "PacBioErrfun.rds", 
#                   "loessErrfun_mod0.rds")){
#   test_by_group <- sequel_unibe_distances |>
#     mutate(dataset.prop = factor(format(dataset.prop, scientific = TRUE), levels = c("1e+04", "1e+05", "1e+06", "1e+07"))) |>
#     group_by(dataset.prop) |>
#     summarise(
#       n         = n(),
#       mean      = mean(.data[[err.func]], na.rm = TRUE),
#       shapiro.p = tryCatch(shapiro.test(.data[[err.func]])$p.value,
#                            error = function(e) NA_real_),
#       p.value   = if_else(shapiro.p > 0.05, tryCatch(t.test(.data[[err.func]], mu = 0, alternative = "greater")$p.value,
#                            error = function(e) NA_real_),
#       tryCatch(wilcox.test(.data[[err.func]], mu = 0, alternative = "greater")$p.value,
#                error = function(e) NA_real_)),
#       .groups = "drop"
#     )|>
#     mutate(err.func = err.func)
  
#   results[[err.func]] <- test_by_group
# }

# all_results <- bind_rows(results) |>
#   mutate(p.adj.bonferroni = p.adjust(p.value, method = "bonferroni"),
#          p.adj.benjamin.hochberg = p.adjust(p.value, method = "BH"))
# all_results <- all_results |> drop_na(shapiro.p)

# all_results |>
#   dplyr::filter(p.adj.bonferroni > 0.01 | p.adj.benjamin.hochberg >0.01) |>
#   print()

# write.table(all_results, file.path(results.path, "nbases_sequel_unibe_summary_test.txt"), row.names = F)

# ##Combination and visualization 
# revio  <- read.table(file.path(results.path, "nbases_revio_unibe_summary_test.txt"),  header = TRUE) |> mutate(platform = "Revio UniBe")
# sequel <- read.table(file.path(results.path, "nbases_sequel_unibe_summary_test.txt"), header = TRUE) |> mutate(platform = "Sequel UniBe")

# all_platforms <- bind_rows(revio, sequel) |>
#   mutate(
#     dataset.prop = factor(format(dataset.prop, scientific = TRUE), levels = c("1e+04", "1e+05", "1e+06", "1e+07")),
#     signif = case_when(
#       p.adj.bonferroni < 0.001 ~ "< 0.001",
#       p.adj.bonferroni < 0.01  ~ "< 0.01",
#       p.adj.bonferroni < 0.05  ~ "< 0.05",
#       TRUE ~ "ns"
#     )
#   )

# ggplot(all_platforms, aes(x = dataset.prop, y = mean, group = err.func, colour = err.func)) +
#   geom_line(position = position_dodge(width = 0.3)) +
#   geom_point(aes(shape = signif), size = 3, position = position_dodge(width = 0.3)) +
#   facet_wrap(~platform, scales = "free_y")  +
#   labs(x = "Dataset proportion (%)", y = "Average distance to the full dataset",
#        color = "Error Function", shape="Adjusted p.value") +
#   theme_bw()
# ggsave(file.path(results.path, "nbases_summary_statistics_test.pdf"))

#Sequence table
seq_tables <- list()
meta.list <- list()
N <- length(names(dada2_results_data))
n=1
for (ref_name in names(dada2_results_data)){
  cat("[",n,"/",N,"] : ", ref_name,"\n")
  seq.tab <- makeSequenceTable(dada2_results_data[[ref_name]])
  info <- dada2_meta[ref_name, ]
  sample.names <- paste0(sub(".rds","",ref_name, fixed = TRUE), "_", rownames(seq.tab))
  meta.list[[ref_name]] <- data.frame(
    fastq.files      = rownames(seq.tab),
    platform         = info$dataset,
    dataset.prop     = info$dataset.prop,
    subset.replicate = info$used.seed,
    error.func       = info$errors.function,
    sample.out       = sample.names
  )
  rownames(seq.tab) <- sample.names
  seq_tables[[ref_name]] <- seq.tab
  n<-n+1
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