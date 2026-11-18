#!/usr/bin/R

#Load libraries
library(phyloseq, verbose = FALSE, quietly = TRUE); packageVersion("phyloseq")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr")
library(cowplot, verbose = FALSE, quietly = TRUE); packageVersion("cowplot")

#Define functions
## Group genera: top N (over the whole phyloseq object) + "Other"; named by taxon
getGenusGroups <- function(ps, rank = "Genus", top.genera = 15) {
  genus <- as.character(tax_table(ps)[, rank])
  genus[is.na(genus) | genus == ""] <- "Unassigned"
  tot <- tapply(taxa_sums(ps), genus, sum)
  keep <- head(names(sort(tot, decreasing = TRUE)), top.genera)
  grp <- ifelse(genus %in% keep, genus, "Other")
  names(grp) <- taxa_names(ps)
  list(groups = grp, levels = c(keep, if (any(grp == "Other")) "Other", "Absent"))
}

## Confusion data of ONE dataset: each subsample vs the full dataset, per error function
## Cell = number of ASVs (matched by sequence) with genus in full (rows) vs genus in subsample (cols)
getConfusion <- function(ps, genus.grp, full.prop) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  meta <- data.frame(sample_data(ps))
  grp <- genus.grp[colnames(otu)]
  
  out <- list()
  for (e in unique(as.character(meta$error.func))) {
    ref.idx <- which(meta$error.func == e & as.character(meta$dataset.prop) == full.prop)
    if (length(ref.idx) == 0) {
      message("  no full dataset for error function '", e, "': skipped")
      next
    }
    ref <- otu[ref.idx[1], ] > 0
    subs <- which(meta$error.func == e & as.character(meta$dataset.prop) != full.prop)
    for (j in subs) {
      sub <- otu[j, ] > 0
      used <- ref | sub
      if (!any(used)) next
      out[[length(out) + 1]] <- data.frame(
        error.func   = e,
        dataset.prop = as.character(meta$dataset.prop[j]),
        full = ifelse(ref[used], grp[used], "Absent"),
        sub  = ifelse(sub[used], grp[used], "Absent")
      ) |> count(error.func, dataset.prop, full, sub, name = "n")
    }
  }
  if (length(out) == 0) return(NULL)
  bind_rows(out) |>
    group_by(error.func, dataset.prop, full, sub) |>   # sums replicates if any
    summarise(n = sum(n), .groups = "drop")
}

## One confusion-matrix figure for ONE error function (rows = dataset, cols = subsample level)
## Both axes and cell values are hidden: only the colour scale encodes the counts
plotConfusion <- function(df, genus.levels, prop.levels) {
  df$full <- factor(df$full, levels = rev(genus.levels))
  df$sub  <- factor(df$sub,  levels = genus.levels)
  df$dataset.prop <- factor(df$dataset.prop, levels = prop.levels)
  
  ggplot(df, aes(x = sub, y = full, fill = n)) +
    geom_tile(colour = "grey80") +
    facet_grid(rows = vars(dataset.label), cols = vars(dataset.prop),
               labeller = labeller(dataset.label = label_wrap_gen(width = 25)),
               drop = FALSE) +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    scale_fill_viridis_c(trans = "log10", name = "# ASVs", direction = -1) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_text(angle = 0),
          panel.grid   = element_blank())
}

#Dataset labels
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

## Automatic label when no match exists in the summary table: cleaned platform name
autoLabel <- function(ds) trimws(gsub("[._]+", " ", ds))

datasets_details <- read.csv("summary/Dataset_Summary-Before_Normalization.csv")
datasets_details <- datasets_details |>
  rowwise() |>
  mutate(
    optional = list({
      o <- c(Library.Type, Segmentation, Chemistry)
      o[!is.na(o) & o != ""]
    }),
    loc.lab = ifelse(is.na(Location)  | Location  == "", "Unknown location",  Location),
    seq.lab = ifelse(is.na(Sequencer) | Sequencer == "", "Unknown sequencer", Sequencer),
    label   = paste0(loc.lab, " (", paste(c(seq.lab, optional), collapse = "/"), ")"),
    key     = norm(paste0(Location, paste(c(Sequencer, optional), collapse = "")))
  ) |>
  ungroup()

#Configuration of the two approaches
approaches <- list(
  prop = list(
    file.prefix = "",
    prop.levels = c("1", "5", "10", "25", "50", "100"),
    full.prop   = "100",
    rep.seed    = 1,
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

top.genera <- 15
tax.rank <- "Genus"
results.path <- "summary"
p <- list()   # p[[approach]][[error.func]]

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
  
  # Same genus grouping for every dataset so that all panels share the same axes
  gg <- getGenusGroups(ps, rank = tax.rank, top.genera = top.genera)
  
  # --- Confusion data of each dataset (subsample vs full) ---
  df.list <- list()
  used.labels <- character(0)
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
    
    # Dataset label (automatic fallback if no match), kept unique across datasets
    idx <- match(norm(ds), datasets_details$key)
    ds.label <- if(is.na(idx)) autoLabel(ds) else datasets_details$label[idx]
    if(is.na(idx)) message("  no label found for '", ds, "': using '", ds.label, "'")
    if(ds.label %in% used.labels) ds.label <- paste0(ds.label, " [", ds, "]")
    used.labels <- c(used.labels, ds.label)
    
    df <- getConfusion(ps.ds, gg$groups, full.prop = cfg$full.prop)
    if(is.null(df)) next
    df$dataset.label <- ds.label
    df.list[[ds]] <- df
  }
  
  if(length(df.list) == 0) next
  df.all <- bind_rows(df.list)
  df.all$dataset.label <- factor(df.all$dataset.label, levels = unique(df.all$dataset.label))
  prop.levels.sub <- setdiff(cfg$prop.levels, cfg$full.prop)
  
  # --- One figure per error function: rows = dataset, cols = subsample level ---
  p[[approach]] <- list()
  message("Working with: ", approach)
  n.ds <- nlevels(df.all$dataset.label)
  for(e in unique(df.all$error.func)){
    df.e <- df.all |> filter(error.func == e)
    p[[approach]][[e]] <- plotConfusion(df.e, gg$levels, prop.levels.sub) +
      ggtitle(paste0("    Error function: ", e))
    ggsave(file.path(results.path,
                     paste0("confusion_", cfg$file.prefix, e, "_all_datasets_plot.pdf")),
           plot = p[[approach]][[e]],
           width = 2.5 * length(prop.levels.sub) + 3, height = 1.5 * n.ds,
           limitsize = FALSE)
  }
}

panel_plot <- list()
for(approach in names(p)){
  panel_plot[[approach]] <- plot_grid(
    plotlist = p[[approach]], 
    ncol = 2, 
    labels = "AUTO", # Automatically labels panels A, B, C, D
    align = "hv"     # Aligns plots horizontally and vertically
  ) 
}

panel_plot[['prop']]

panel_plot[['']]

message("Done.")