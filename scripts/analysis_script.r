#!/usr/bin/R

#Load libraries
library(dada2); packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")
library(phyloseq); packageVersion("phyloseq")

#Check whether the number of parameters is 2
if(length(args)!=2){
  print("This script takes the following parameters:")
  print("   (1) Path to input directory")
  print("   (2) Path to output directory")
  stop("Error: Required arguments are not provided.", call.=FALSE)
}

#Setup inputs and outputs 
path.input <- args[1]
path.output <- args[2]
path.rds <- file.path(path.output, "RDS")
path.figures <- file.path(path.output, "Figure")
fastq.input <- list.files(path.input, pattern="fastq.gz", full.names=TRUE)

#Quality Control
plotQualityProfile(fastq.input)
lens.fn <- lapply(fastq.input, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)