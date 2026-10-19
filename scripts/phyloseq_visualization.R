#!/usr/bin/R

#Load libraries
library(phyloseq, verbose = FALSE, quietly = TRUE); packageVersion("phyloseq")
library(Biostrings, verbose = FALSE, quietly = TRUE); packageVersion("Biostrings")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr")

#Define functions
## Compute the relative abundances of the top taxa of ONE dataset (long format)
getTaxaDistrib <- function(ps, top = 20, fill = "Family") {
  top.taxa <- head(names(sort(taxa_sums(ps), decreasing = TRUE)), top)
  ps.freq <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
  ps.freq <- prune_taxa(top.taxa, ps.freq)
  
  psmelt(ps.freq) %>%
    group_by(dataset.prop, error.func, .data[[fill]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
}

## Plot all datasets in one figure: rows = dataset, columns = error function
plotTaxaDistrib <- function(df, fill = "Family", xlab = "Dataset proportion (%)") {
  ggplot(df, aes(x = factor(dataset.prop), y = Abundance, fill = .data[[fill]])) +
    geom_col(position = "stack", colour = "black") +
    facet_grid(rows = vars(dataset.label), cols = vars(error.func),
               labeller = labeller(dataset.label = label_wrap_gen(width = 25))) +
    labs(x = xlab, y = "Abondance", fill = fill) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 0))
}

#Dataset labels
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

datasets_details <- read.csv("summary/Dataset_Summary-Before_Normalization.csv")
datasets_details <- datasets_details |>
  rowwise() |>
  mutate(
    optional = list({
      o <- c(Library.Type, Segmentation, Chemistry)
      o[!is.na(o) & o != ""]
    }),
    label = paste0(Location, " (", paste(c(Sequencer, optional), collapse = "/"), ")"),
    key   = norm(paste0(Location, paste(c(Sequencer, optional), collapse = "")))
  ) |>
  ungroup()

## Column of datasets_details that matches the `platform` column of the metadata
## (ADAPT if needed)
details.key <- "Platform"

#Configuration of the two approaches
approaches <- list(
  prop = list(
    file.prefix = "",
    prop.levels = c("1", "5", "10", "25", "50", "100"),
    full.prop   = "100",
    rep.seed    = 41,
    xlab        = "Dataset proportion (%)",
    input_data  = "summary/proportion"
  ),
  nbases = list(
    file.prefix = "nbases_",
    prop.levels = c("1e+04", "1e+05", "1e+06", "1e+07", "1e+08"),
    full.prop   = "1e+08",
    rep.seed    = NULL,
    xlab        = "Number of bases for learning errors",
    input_data  = "summary/nbases"
  )
)

top <- 100
results.path <- "summary"
p <- list()
i<-1
for(approach in names(approaches)){
  cfg <- approaches[[approach]]
  path.seqtab <- file.path(cfg$input_data, paste0(cfg$file.prefix, "all_dataset_sequence_table_nochim.rds"))
  path.taxa   <- file.path(cfg$input_data, paste0(cfg$file.prefix, "all_dataset_taxonomy_assignment.rds"))
  path.meta   <- file.path(cfg$input_data, paste0(cfg$file.prefix, "all_metadata.csv"))
  if(!all(file.exists(path.seqtab, path.taxa, path.meta))){
    message("[", approach, "] input files not found in ", cfg$input_data, ": skipped")
    next
  }
  message("[", approach, "] loading data")
  
  combined.seq_tables.nochim <- readRDS(path.seqtab)
  taxaAssign <- readRDS(path.taxa)
  all.meta.data <- read.csv(path.meta, colClasses = c(dataset.prop = "character"))
  rownames(all.meta.data) <- all.meta.data[,1]
  all.meta.data[,1] <- NULL
  all.meta.data$dataset.prop[is.na(all.meta.data$dataset.prop)] <- cfg$full.prop
  all.meta.data$dataset.prop <- factor(all.meta.data$dataset.prop, levels = cfg$prop.levels)
  
  ps <- phyloseq(otu_table(combined.seq_tables.nochim, taxa_are_rows = FALSE),
                 sample_data(all.meta.data),
                 tax_table(taxaAssign))
  
  # --- Collect the distribution of each dataset ---
  df.list <- list()
  for(ds in unique(all.meta.data$platform)){
    message("[", approach, "] ", ds)
    
    meta.ds <- data.frame(sample_data(ps))
    keep <- meta.ds$platform == ds
    if(!is.null(cfg$rep.seed)){
      keep <- keep & (is.na(meta.ds$subset.replicate) | meta.ds$subset.replicate == cfg$rep.seed)
      missing.prop <- setdiff(unique(as.character(meta.ds$dataset.prop[meta.ds$platform == ds])),
                              unique(as.character(meta.ds$dataset.prop[keep])))
      if(length(missing.prop) > 0){
        message("  seed ", cfg$rep.seed, " not available for proportion(s): ", paste(missing.prop, collapse = ", "))
      }
    }
    if(!any(keep)){
      message("  no samples left: skipped")
      next
    }
    ps.ds <- prune_samples(keep, ps)
    ps.ds <- prune_taxa(taxa_sums(ps.ds) > 0, ps.ds)
    
    dna <- Biostrings::DNAStringSet(taxa_names(ps.ds))
    names(dna) <- taxa_names(ps.ds)
    ps.ds <- merge_phyloseq(ps.ds, dna)
    taxa_names(ps.ds) <- paste0("ASV", seq(ntaxa(ps.ds)))
    
    # Dataset label (falls back to the platform name if no match)
    idx <- match(norm(ds), datasets_details$key)
    ds.label <- if(is.na(idx)) ds else datasets_details$label[idx]
    if(is.na(idx)) message("  no label found for '", ds, "'")
    
    df <- getTaxaDistrib(ps.ds, top = top, fill = "Species")
    df$dataset.label <- ds.label
    df.list[[ds]] <- df
  }
  
  if(length(df.list) == 0) next
  df.all <- bind_rows(df.list)
  df.all$dataset.label <- factor(df.all$dataset.label, levels = unique(df.all$dataset.label))
  
  # --- One figure: rows = dataset, columns = error function ---

  p[[i]] <- plotTaxaDistrib(df.all, fill = "Species", xlab = cfg$xlab)
  n.ds <- nlevels(df.all$dataset.label)
  n.err <- length(unique(df.all$error.func))
  ggsave(file.path(results.path, paste0("taxonomic_distrib_top_", top, "_", cfg$file.prefix, "all_datasets_plot.pdf")),
         plot = p[[i]], width = 3 * n.err + 4, height = 3 * n.ds, limitsize = FALSE)
  i<-i+1
}

message("Done.")
quit(save = "no", status = 0)