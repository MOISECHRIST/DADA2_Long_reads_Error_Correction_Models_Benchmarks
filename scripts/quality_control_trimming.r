#!/usr/bin/R

#Load libraries
library(dada2, verbose = FALSE, quietly = TRUE); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")

#Check whether the number of parameters is 2
args <- commandArgs(trailingOnly=TRUE)

if(length(args)!=2){
  print("This script takes the following parameters:")
  print("   (1) Path to input directory")
  print("   (2) Path to output directory")
  stop("Error: Required arguments are not provided.", call.=FALSE)
}

#Setup inputs and outputs 
path.input <- args[1]
path.output <- args[2]
path.figures <- file.path(path.output, "Figure")
fastq.input <- list.files(path.input, pattern="fastq.gz", full.names=TRUE)
dir.create(path.figures, recursive = TRUE, showWarnings = FALSE)

#Quality Control
pdf(file.path(path.figures, "quality_profile.pdf"))
print(plotQualityProfile(fastq.input))
dev.off()

#Reads length distribution
lens.fn <- lapply(fastq.input, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
pdf(file.path(path.figures, "hist_len_plot.pdf"))
hist(lens, 100)
dev.off()


#Filter en trimming
filts <- file.path(path.output, "Filtered", basename(fastq.input))
track <- filterAndTrim(fastq.input, filts, minQ=3, minLen=1000, maxLen=1600
  maxN=0, rm.phix=FALSE, maxEE=2, multithread = TRUE)
track <- as.data.frame(track)
track$sample <- rownames(track)
rownames(track) <- NULL 

ggplot(data=track)+
  geom_col(aes(x=reads.in, y=sample, fill="Reads in"), position="dodge")+
  geom_col(aes(x=reads.out, y=sample, fill="Reads out"), position="dodge")+
  scale_fill_manual(name="Legend", values=c("Reads in"="darkblue", "Reads out"="darkgreen"))+
  labs(x="Number of reads")
ggsave(filename = file.path(path.figures, "track_filter-trim_plot.pdf"))
write.csv(track, file.path(path.figures, "track_filter-trim.csv"))


























