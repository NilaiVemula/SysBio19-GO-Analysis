#!/usr/bin/env Rscript

require(optparse)

setwd("../logs")
print(getwd())

source("../scripts/04-functional_analysis/FunctionalEnrichmentAnalysis.R")

## parse command line args
parseCommandLineArgs <- function() {
    options <- list(
                    make_option(c("-r", "--runName"),
                                help="run name",
                                type = "character")
                    ) # end optionList definition

    parse_args(OptionParser(option_list = options))
}

args <- parseCommandLineArgs()
run_name = args$runName

dir <- getwd()
file <- paste("../results/post-processing/",run_name,"/fa_annotation/gene_expression_and_module_membership.txt",sep="")

gene_expression_and_module_membership <- read.delim(file)
dt <- data.frame(gene_expression_and_module_membership$Gene,gene_expression_and_module_membership$Module,gene_expression_and_module_membership$Metamodule)
colnames(dt) <- c("ensemblid","module_name","metamodule")


output_folder <- paste("../results/post-processing/", run_name, "/functional-enrichment/", sep="")
dir.create(file.path(output_folder), showWarnings = FALSE)

module_list <- unique(dt$module_name)

for(key_module
    in module_list){

module_genes_df   <- dt[dt$module_name == key_module, ]

f <- geneEnrichmentAnalysis(
    output_folder = output_folder,
    gene_list = module_genes_df,
    title = paste0("module_", key_module)
)
}


metamodule_list <- unique(dt$metamodule)

for(key_metamodule
    in metamodule_list){

  metamodule_genes_df <- dt[dt$metamodule == key_metamodule, ]
  f <- geneEnrichmentAnalysis(
      output_folder = output_folder,
      gene_list = metamodule_genes_df,
      title = paste0("metamodule_", key_metamodule)
  )

}
