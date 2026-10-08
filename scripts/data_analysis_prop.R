#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(tidyr, verbose = FALSE, quietly = TRUE); packageVersion("tidyr") 
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr") 

#Constants
## Proportion (%) of the dataset used for learning errors. "100" = full dataset, 
## i.e. runs whose directory name has no "_<prop>_<seed>" suffix.
prop.levels <- c("1", "5", "10", "25", "50", "100")
full.prop <- "100"
## Resource metrics written by track_resources() in dada2_analysis.R
resource.cols <- c("wall_sec", "cpu_sec", "avg_cores_used", "threads_available",
                   "cpus_allocated", "peak_ram_gb", "end_ram_gb")

#Define functions

## Parse <dataset>[_<prop>_<seed>]/RDS/<prefix>_<errors.function>.<ext>
## Dataset names now contain underscores (e.g. LIB_16S_KINNEX_SEGMENTED_REVIO_SPRQ_NX), 
## so we cannot split on "_" anymore: the proportion and seed suffix is matched 
## from the end of the run directory name.
parseRunPath <- function(file_path){
  run  <- basename(dirname(dirname(file_path)))
  file <- basename(file_path)
  pat  <- "^(.+?)_(1|5|10|25|50)_([0-9]+)$"
  has.suffix <- grepl(pat, run, perl = TRUE)
  
  errors.function <- sub("^(dada_results|learn_error_data|execution_time)_", "", file)
  errors.function <- sub("\\.(rds|csv)$", "", errors.function)
  
  data.frame(
    ref_name        = paste0(run, "_", file),
    run             = run,
    dataset         = ifelse(has.suffix, sub(pat, "\\1", run, perl = TRUE), run),
    dataset.prop    = ifelse(has.suffix, sub(pat, "\\2", run, perl = TRUE), full.prop),
    used.seed       = ifelse(has.suffix, sub(pat, "\\3", run, perl = TRUE), NA_character_),
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

readMultiCSV <- function(list_of_paths){
  res <- list()
  for(file_path in list_of_paths){
    info <- parseRunPath(file_path)
    df <- read.csv(file_path)
    df$X <- NULL
    #Drop files where learnError is the only process logged (failed runs)
    if(all(df$process_name == "learnError")) next
    df$platform        <- info$dataset
    df$dataset.prop    <- info$dataset.prop
    df$used.seed       <- info$used.seed
    df$errors.function <- info$errors.function
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
  return(
    res
  )
}

## Distance between the error matrix learned on the full dataset and the one 
## learned on a subset, for every error function available in the dataset.
## Returns a wide data.frame (one column "<errors.function>.rds" per function).
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
      dataset.prop = factor(dataset.prop, levels = prop.levels)
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

## Generic resource plot: replicates as dodged points, median as line, 
## all datasets in the same panel (one panel per error function)
plotResourcePlot.prop <- function(execution.time.long, metric, ylab){
  execution.time.long |>
    ggplot()+
    geom_boxplot(aes(x=dataset.prop, y=.data[[metric]], colour = process_name), alpha=0.6)+
    facet_grid(cols=vars(errors.function), rows = vars(platform), scales = "free",
               labeller = labeller(errors.function = label_wrap_gen(width = 100),
                                   platform = label_wrap_gen(width = 20)))+
    labs(x="Dataset proportion (%)",
         y=ylab,
         colour = "Process Name") + theme_bw()
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

## Is the distance to the full dataset significantly > 0 at each proportion?
## (t-test if Shapiro-Wilk does not reject normality, Wilcoxon otherwise).
## Needs replicates (seeds); returns NULL if there is nothing to test.
runDistanceTests <- function(distances){
  results <- list()
  for(err.func in grep("\\.rds$", names(distances), value = TRUE)){
    test_by_group <- distances |>
      #The full dataset is at distance 0 from itself: nothing to test
      filter(dataset.prop != as.numeric(full.prop)) |>
      mutate(dataset.prop = factor(dataset.prop, levels = setdiff(prop.levels, full.prop))) |>
      group_by(dataset.prop) |>
      summarise(
        n         = n(),
        mean      = mean(.data[[err.func]], na.rm = TRUE),
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
path.alldataset <- file.path(list.files("results_prop", full.names = T), "RDS")
paths_learn_errors_data <- list.files(path.alldataset, pattern = "learn_error", full.names = T) 
paths_dada2_results_data <- list.files(path.alldataset, pattern = "dada_results", full.names = T) 
paths_exec_times <- list.files(path.alldataset, pattern = "execution_time", full.names = T)
ref.db <- file.path("refSeq/SILVA-v138.2-16s/silva_nr99_v138.2_toSpecies_trainset.fa.gz")

results.path <- "summary"
dir.create(results.path, showWarnings = F)

#Load data
learn_errors_data <- readMultiRDS(paths_learn_errors_data)
dada2_results_data <- readMultiRDS(paths_dada2_results_data)
exec_times_data <- readMultiCSV(paths_exec_times)

#Run metadata (dataset / proportion / seed / error function), indexed by ref_name
dada2_meta <- bind_rows(lapply(paths_dada2_results_data, parseRunPath)) |> as.data.frame()
rownames(dada2_meta) <- dada2_meta$ref_name
datasets <- unique(dada2_meta$dataset)
message(length(datasets), " datasets found:\n  ", paste(datasets, collapse = "\n  "))

#Execution time and resource usage
## Execution time
p <- plotExecutionTimePlot.prop(exec_times_data)
saveFacetPlot(p, "Execution_Time_boxplot.pdf", exec_times_data, by.platform = FALSE)

## Total CPU time
p <- plotResourcePlot.prop(exec_times_data, "cpu_sec", "CPU Time (s)")
saveFacetPlot(p, "CPU_Time_plot.pdf", exec_times_data, by.platform = FALSE)

## Average number of busy cores (cpu_sec / wall_sec)
p <- plotResourcePlot.prop(exec_times_data, "avg_cores_used", "Average cores used")
saveFacetPlot(p, "Cores_Used_plot.pdf", exec_times_data, by.platform = FALSE)

## CPU efficiency: average busy cores / cores available
p <- plotResourcePlot.prop(exec_times_data, "cpu_efficiency", "CPU efficiency (avg cores used / cores available)")
saveFacetPlot(p, "CPU_Efficiency_plot.pdf", exec_times_data, by.platform = FALSE)

## Peak RAM
## NB: peak is per step only if /proc/self/clear_refs was writable; otherwise 
## VmHWM is cumulative and Denoising peak >= learnError peak.
p <- plotResourcePlot.prop(exec_times_data, "peak_ram_gb", "Peak RAM (GB)")
saveFacetPlot(p, "Peak_RAM_plot.pdf", exec_times_data, by.platform = FALSE)

## Summary table (mean and sd over replicates/seeds)
resource_summary <- exec_times_data |>
  group_by(platform, errors.function, dataset.prop, process_name) |>
  summarise(
    n_runs = n(),
    across(c(duration, cpu_sec, avg_cores_used, cpu_efficiency, peak_ram_gb, end_ram_gb),
           list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))),
    .groups = "drop"
  )
write.csv(resource_summary, file.path(results.path, "resource_summary.csv"), row.names = F)

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
  ggsave(file.path(results.path, paste0(ds, "_dataset_distance_plot.pdf")), plot = p)
  p <- plotDistanceBoxplot.prop(ds.distances.long)
  ggsave(file.path(results.path, paste0(ds, "_dataset_distance_boxplot.pdf")), plot = p)
}

#Combination of datasets
combined_distances.long <- bind_rows(distances_long_by_dataset)
write.csv(combined_distances.long, file.path(results.path, "all_distances.csv"), row.names = F)

p <- combined_distances.long |>
  ggplot()+
  geom_boxplot(aes(x=dataset.prop, y=distances), outliers = F)+
  geom_jitter(aes(x=dataset.prop, y=distances, colour = dataset.prop), alpha=0.6)+
  facet_grid(cols=vars(errors.function), rows = vars(platform), scales = "free",
             labeller = labeller(errors.function = label_wrap_gen(width = 100),
                                 platform = label_wrap_gen(width = 20)))+
  labs(x="Dataset proportion (%)",
       y="Distance to the full dataset") + 
  labs(colour="Dataset\nproportion (%)") + theme_bw()
saveFacetPlot(p, "Combined_dataset_distance_boxplot.png", combined_distances.long)

dodge <- position_dodge(width = 0.3)
p <- combined_distances.long |>
  ggplot(aes(x=dataset.prop, y=distances, colour = platform, shape = platform, group = platform))+
  stat_summary(fun = median, geom = "line", position = dodge, size=1)+
  #geom_point(position = dodge, alpha=0.4)+
  scale_shape_manual(values = rep_len(c(16, 17, 15, 3, 7, 8, 4, 0, 1, 2, 5, 6, 9, 10, 11, 12, 13, 14),
                                      n_distinct(combined_distances.long$platform)))+
  facet_grid(cols=vars(errors.function), scales = "free",
             labeller = labeller(errors.function = label_wrap_gen(width = 100)))+
  labs(x="Dataset proportion (%)",
       y="Distance to the full dataset",
       colour="Dataset", shape="Dataset")+theme_bw()
saveFacetPlot(p, "Combined_dataset_distance_plot.pdf", combined_distances.long, by.platform = FALSE)


#Statistics tests
## One test table per dataset (written as <dataset>_summary_test.txt, lower case)
test_results <- list()
for(ds in names(distances_by_dataset)){
  ds.results <- runDistanceTests(distances_by_dataset[[ds]])
  if(is.null(ds.results) || nrow(ds.results) == 0){
    message("[", ds, "] not enough replicates for statistical tests: skipped")
    next
  }
  write.table(ds.results, file.path(results.path, paste0(tolower(ds), "_summary_test.txt")), row.names = F)
  ds.results |>
    dplyr::filter(p.adj.bonferroni > 0.01 | p.adj.benjamin.hochberg >0.01) |>
    print()
  test_results[[ds]] <- ds.results
}

##Combination and visualization 
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

if(nrow(all_platforms) > 0){
  p <- ggplot(all_platforms, aes(x = dataset.prop, y = mean, group = err.func, colour = err.func)) +
    geom_line(position = position_dodge(width = 0.3)) +
    geom_point(aes(shape = signif), size = 3, position = position_dodge(width = 0.3)) +
    facet_wrap(~platform, scales = "free_y", ncol = 3,
               labeller = labeller(platform = label_wrap_gen(width = 30)))  +
    labs(x = "Dataset proportion (%)", y = "Average distance to the full dataset",
         color = "Error Function", shape="Adjusted p.value") +
    theme_bw()
  ggsave(file.path(results.path, "summary_statistics_test.png"), plot = p,
         width = 14, height = 4 * ceiling(n_distinct(all_platforms$platform) / 3), limitsize = FALSE)
}

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
write.csv(all.meta.data,file.path(results.path, "all_metadata.csv"))
combined.seq_tables <- mergeSequenceTables(tables=seq_tables)
saveRDS(combined.seq_tables,file.path(results.path, "all_dataset_sequence_table.rds"))

#Assign taxonomy
combined.seq_tables.nochim <- removeBimeraDenovo(combined.seq_tables, method="consensus", multithread=TRUE)
saveRDS(combined.seq_tables.nochim,file.path(results.path, "all_dataset_sequence_table_nochim.rds"))
taxaAssign <- assignTaxonomy(combined.seq_tables.nochim, ref.db,
                             multithread = T)
saveRDS(taxaAssign,file.path(results.path, "all_dataset_taxonomy_assignment.rds"))


message("Done.")
quit(save = "no", status = 0)