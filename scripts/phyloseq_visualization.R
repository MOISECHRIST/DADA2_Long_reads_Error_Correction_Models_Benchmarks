#!/usr/bin/R

#Load libraries
library(phyloseq, verbose = FALSE, quietly = TRUE); packageVersion("phyloseq")
library(Biostrings, verbose = FALSE, quietly = TRUE); packageVersion("Biostrings")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")
library(dplyr, verbose = FALSE, quietly = TRUE); packageVersion("dplyr") 

#Define functions
plotTaxaDistrib <- function(ps, top = 20, fill = "Family", xlab = "Dataset proportion (%)") {
  top.taxa <- head(names(sort(taxa_sums(ps), decreasing = TRUE)), top)
  ps.freq <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
  ps.freq <- prune_taxa(top.taxa, ps.freq)
  
  df <- psmelt(ps.freq)
  
  df_grouped <- df %>%
    group_by(dataset.prop, error.func, .data[[fill]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  ggplot(df_grouped, aes(x = factor(dataset.prop), y = Abundance, fill = .data[[fill]])) +
    geom_col(position = "stack", colour="black") +
    facet_grid(cols = vars(error.func)) +
    labs(x = xlab, y = "Abondance", fill = fill) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#Configuration of the two approaches
## prop   : subsets of the dataset (results_prop, produced by data_analysis_prop.R)
##          runs <dataset>_<prop>_<seed>; one representative seed is kept per proportion
## nbases : subsets by number of bases (results_nbases, produced by data_analysis_nbases.R)
##          one run per number of bases, so no seed to choose
## An approach whose files are not found in `summary/` is skipped.
approaches <- list(
  prop = list(
    file.prefix = "",
    prop.levels = c("1", "5", "10", "25", "50", "100"),
    full.prop   = "100",
    rep.seed    = 41,     # seed kept as the representative replicate (NULL = keep all)
    xlab        = "Dataset proportion (%)"
  ),
  nbases = list(
    file.prefix = "nbases_",
    prop.levels = c("1e+04", "1e+05", "1e+06", "1e+07", "1e+08"),
    full.prop   = "1e+08",
    rep.seed    = NULL,
    xlab        = "Number of bases for learning errors"
  )
)

top <- 100

#From here I start following the tutorial: https://benjjneb.github.io/dada2/tutorial.html 
#On Bonus: Handoff to phyloseq

results.path <- "summary"

for(approach in names(approaches)){
  cfg <- approaches[[approach]]
  path.seqtab <- file.path(results.path, paste0(cfg$file.prefix, "all_dataset_sequence_table_nochim.rds"))
  path.taxa   <- file.path(results.path, paste0(cfg$file.prefix, "all_dataset_taxonomy_assignment.rds"))
  path.meta   <- file.path(results.path, paste0(cfg$file.prefix, "all_metadata.csv"))
  if(!all(file.exists(path.seqtab, path.taxa, path.meta))){
    message("[", approach, "] input files not found in ", results.path, ": skipped")
    next
  }
  message("[", approach, "] loading data")
  
  combined.seq_tables.nochim <- readRDS(path.seqtab)
  taxaAssign <- readRDS(path.taxa)
  #dataset.prop must stay a character ("1e+04" would otherwise be read as a number)
  all.meta.data <- read.csv(path.meta, colClasses = c(dataset.prop = "character"))
  rownames(all.meta.data) <- all.meta.data[,1]
  all.meta.data[,1] <- NULL
  all.meta.data$dataset.prop[is.na(all.meta.data$dataset.prop)] <- cfg$full.prop
  all.meta.data$dataset.prop <- factor(all.meta.data$dataset.prop, levels = cfg$prop.levels)
  
  ps <- phyloseq(otu_table(combined.seq_tables.nochim, taxa_are_rows=FALSE),
                 sample_data(all.meta.data),
                 tax_table(taxaAssign))
  
  # --- One analysis per dataset ---
  for(ds in unique(all.meta.data$platform)){
    ds.label <- paste0(cfg$file.prefix, tolower(ds))
    message("[", approach, "] ", ds)
    
    #Keep this dataset only and, for prop, one representative seed per proportion 
    #(the full dataset has no seed)
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
    #Remove the ASVs that belong to the other datasets
    ps.ds <- prune_taxa(taxa_sums(ps.ds) > 0, ps.ds)
    
    dna <- Biostrings::DNAStringSet(taxa_names(ps.ds))
    names(dna) <- taxa_names(ps.ds)
    ps.ds <- merge_phyloseq(ps.ds, dna)
    taxa_names(ps.ds) <- paste0("ASV", seq(ntaxa(ps.ds)))
    
    p <- plot_richness(ps.ds, x="dataset.prop", measures=c("Shannon", "Simpson"), color="error.func")+
      labs(x=cfg$xlab,
           color="Error Function")
    ggsave(file.path(results.path, paste0("alpha_divsersity_", ds.label, "_plot.pdf")), plot = p)
    # Warning message:
    #   In estimate_richness(physeq, split = TRUE, measures = measures) :
    #   The data you have provided does not have
    # any singletons. This is highly suspicious. Results of richness
    # estimates (for example) are probably unreliable, or wrong, if you have already
    # trimmed low-abundance taxa from the data.
    # 
    # We recommended that you find the un-trimmed data and retry.
    
    p <- plotTaxaDistrib(ps.ds, top = top, fill = "Species", xlab = cfg$xlab)
    ggsave(file.path(results.path, paste0("taxonomic_distrib_top_", top, "_", ds.label, "_plot.pdf")), plot = p)
  }
}

message("Done.")
quit(save = "no", status = 0)