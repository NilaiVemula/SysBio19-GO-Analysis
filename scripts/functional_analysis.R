#!/usr/bin/env Rscript

library(tidyverse)

source("scripts/FunctionalEnrichmentAnalysis.R")

genes_df <- read_csv("data/gene_list.csv")

output_folder <- "outputs/"

attractor_list <- unique(genes_df$Attractor_State)

for(key_attractor
    in attractor_list){
  
  attractor_genes_df <- genes_df[genes_df$Attractor_State == key_attractor, ]
  
  f <- geneEnrichmentAnalysis(
    output_folder = output_folder,
    gene_list = attractor_genes_df,
    title = paste0("attractor_", key_attractor)
  )
}



