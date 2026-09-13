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
for(prop in c("1", "5", "10", "25", "50")){
  for(seed in 1:10){
    ps <- prune_samples(!startsWith(sample_names(ps), paste0("Revio_UniBe_",prop, "_",seed,"_")), ps)
  }
}

for(prop in c("1", "5", "10", "25", "50")){
  for(seed in 1:10){
    ps <- prune_samples(!startsWith(sample_names(ps), paste0("Sequel_UniBe_",prop, "_",seed,"_")), ps)
  }
}

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

top <- 900
top.taxa <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:top]
ps.freq <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.freq <- prune_taxa(top.taxa, ps.freq)
plot_bar(ps.freq, x="dataset.prop", fill="Family", ) + 
  facet_grid(rows = vars(platform), scales = "free")
ggsave(file.path(results.path,paste0("taxonomic_distribution_top_",top,"_plot.pdf")))
