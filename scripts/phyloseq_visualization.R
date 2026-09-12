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

ps <- phyloseq(otu_table(combined.seq_tables.nochim, taxa_are_rows=FALSE),
               sample_data(all.meta.data),
               tax_table(taxaAssign))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

top.taxa <- names(sort(taxa_sums(ps), decreasing=TRUE))
ps.freq <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.freq <- prune_taxa(top.taxa, ps.freq)
plot_bar(ps.freq, x="Day", fill="Family") + 
  facet_grid(cols=vars(error.func), rows = vars(platform), scales = "free")