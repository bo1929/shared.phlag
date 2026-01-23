if (!require("phylter")) {
  install.packages("phylter",repos = "http://cran.us.r-project.org")
}
if (!require(ape)) {
  install.packages(ape,repos = "http://cran.us.r-project.org")
}
# setwd("./phylter/R")
# source("phylter.R")
library("phylter")
library(ape)
library(tidyr)
library(dplyr)
library(purrr)

rep="01"
lam="0"
trees_dir <- paste0("./aln200t2m-gene-vs-gene/", rep,"/modifiedTrees_100_lam_", lam, ".trees")

trees <- ape::read.tree(trees_dir)
results <- phylter(trees, InitialOnly = FALSE, gene.names = c(1:length(trees)))
num_sp <- row.names(results$Final$F)
num_gn <- names(results$Final$PartialF)

#gene-species outliers
outliers <- as.data.frame(results$Final$Outliers)
names(outliers) <- c("gene", "species")
head(outliers)

#full gene and species outliers
full_species <- results$Final$CompleteOutliers$ComplOutSP

full_genes <- results$Final$CompleteOutliers$ComplOutGN
