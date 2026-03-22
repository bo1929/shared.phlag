#!/usr/bin/env Rscript
# install.packages("tidyr",repos="http://cran.us.r-project.org")
# install.packages("dplyr",repos="http://cran.us.r-project.org")
# install.packages("purrr",repos="http://cran.us.r-project.org")
# install.packages("ape",repos="http://cran.us.r-project.org")
# install.packages("phylter",repos="http://cran.us.r-project.org")
library(dplyr)
library(tidyr)
library(purrr)
library(ape)
library(phylter)
print("Loaded.")

args=commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Usage: phylter.R <k> <genetree_file> <output_file>")

}
print(args)

input_path <- args[2]
output_path <- args[3]
k <- as.numeric(args[1])

trees <- ape::read.tree(input_path)
results <- phylter(trees, parallel=FALSE, verbose=TRUE, test.island=FALSE, k=k, InitialOnly=FALSE, gene.names=c(1:length(trees)))

# num_sp <- row.names(results$Final$F)
# num_gn <- names(results$Final$PartialF)

#gene-species outliers
# outliers <- as.data.frame(results$Final$Outliers)
# names(outliers) <- c("gene", "species")

# full_species <- results$Final$CompleteOutliers$ComplOutSP
# full_genes <- results$Final$CompleteOutliers$ComplOutGN

write.phylter(results, output_path)
