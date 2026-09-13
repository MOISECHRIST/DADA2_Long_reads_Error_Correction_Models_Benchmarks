#!/usr/bin/R

#Load libraries
library(phyloseq, verbose = FALSE, quietly = TRUE); packageVersion("phyloseq")
library(Biostrings, verbose = FALSE, quietly = TRUE); packageVersion("Biostrings")
library(ggplot2, verbose = FALSE, quietly = TRUE); packageVersion("ggplot2")

#From here I start following the tutorial: https://benjjneb.github.io/dada2/tutorial.html 
#On Bonus: Handoff to phyloseq

results.path <- "summary"
combined.seq_tables.nochim <- readRDS(file.path(results.path, "all_dataset_sequence_table_nochim.rds"))
taxaAssign <- readRDS(file.path(results.path, "all_dataset_taxonomy_assignment.rds"))
all.meta.data <- read.csv(file.path(results.path,"all_metadata.csv"))
rownames(all.meta.data) <- all.meta.data[,1]
all.meta.data[,1] <- NULL
all.meta.data$dataset.prop <- factor(all.meta.data$dataset.prop, levels = c("1", "5", "10", "25", "50", "100"))
all.meta.data$dataset.prop[is.na(all.meta.data$dataset.prop)] <- "100"

ps <- phyloseq(otu_table(combined.seq_tables.nochim, taxa_are_rows=FALSE),
               sample_data(all.meta.data),
               tax_table(taxaAssign))

#Remove seed replicates
ps.revio_unibe <- ps
for(prop in c("1", "5", "10", "25", "50")){
  for(seed in 1:10){
    ps.revio_unibe <- prune_samples(!startsWith(sample_names(ps.revio_unibe), paste0("Revio_UniBe_",prop, "_",seed,"_")), ps.revio_unibe)
    ps.revio_unibe <- prune_samples(!startsWith(sample_names(ps.revio_unibe), "Sequel_UniBe_"), ps.revio_unibe)
  }
}

dna <- Biostrings::DNAStringSet(taxa_names(ps.revio_unibe))
names(dna) <- taxa_names(ps.revio_unibe)
ps.revio_unibe <- merge_phyloseq(ps.revio_unibe, dna)
taxa_names(ps.revio_unibe) <- paste0("ASV", seq(ntaxa(ps)))

# --- Revio UniBe---
ps.revio_unibe <- ps
for(prop in c("1", "5", "10", "25", "50")){
  for(seed in 1:10){
    ps.revio_unibe <- prune_samples(!startsWith(sample_names(ps.revio_unibe), paste0("Revio_UniBe_",prop, "_",seed,"_")), ps.revio_unibe)
    ps.revio_unibe <- prune_samples(!startsWith(sample_names(ps.revio_unibe), "Sequel_UniBe_"), ps.revio_unibe)
  }
}
dna <- Biostrings::DNAStringSet(taxa_names(ps.revio_unibe))
names(dna) <- taxa_names(ps.revio_unibe)
ps.revio_unibe <- merge_phyloseq(ps.revio_unibe, dna)
taxa_names(ps.revio_unibe) <- paste0("ASV", seq(ntaxa(ps.revio_unibe)))

# --- Sequel UniBe---
ps.sequel_unibe <- ps
for(prop in c("1", "5", "10", "25", "50")){
  for(seed in 1:10){
    ps.sequel_unibe <- prune_samples(!startsWith(sample_names(ps.sequel_unibe), paste0("Sequel_UniBe_",prop, "_",seed,"_")), ps.sequel_unibe)
    ps.sequel_unibe <- prune_samples(!startsWith(sample_names(ps.sequel_unibe), "Revio_UniBe_"), ps.sequel_unibe)
  }
}
dna <- Biostrings::DNAStringSet(taxa_names(ps.sequel_unibe))
names(dna) <- taxa_names(ps.sequel_unibe)
ps.sequel_unibe <- merge_phyloseq(ps.sequel_unibe, dna)
taxa_names(ps.sequel_unibe) <- paste0("ASV", seq(ntaxa(ps.sequel_unibe)))

plotTaxaDistrib <- function(ps, top = 20, fill="Family"){
  top.taxa <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:top]
  ps.freq <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps.freq <- prune_taxa(top.taxa, ps.freq)
  plot_bar(ps.freq, x="dataset.prop", fill=fill, ) + 
    facet_grid(cols = vars(error.func)) +
    labs(x="Dataset proportion (%)")
}

top <- 30
plotTaxaDistrib(ps.revio_unibe, top = top)
ggsave(file.path(results.path,paste0("taxonomic_distrib_top_",top,"_revio_unibe_plot.pdf")))
plotTaxaDistrib(ps.sequel_unibe, top = top)
ggsave(file.path(results.path,paste0("taxonomic_distrib_top_",top,"_sequel_unibe_plot.pdf")))
