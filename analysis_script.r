#!/usr/bin/R

BiocManager::install("phyloseq", version = "3.16")

library(dada2);packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")
library(phyloseq); packageVersion("phyloseq")
