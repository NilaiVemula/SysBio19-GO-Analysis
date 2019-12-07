library(clusterProfiler)
library("org.Mm.eg.db")


geneEnrichmentAnalysis <- function(output_folder, gene_list,  title){
    
    # fake data
    #output_folder = "results/WT-HFD-females__vs__WT-HFD-males/functional-enrichment/"
    #gene_list = genes_upregulated_df
    #title = "WT-HFD-females__vs__WT-HFD-males__Upregulated"

    
    #  
    # Gene Enrichment Analysis  --------------------------------------------------
    # 
    

    
    # GO Classification
    # 

    # http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#go-analysis
    # groupGO -  Functional Profile of a gene set at specific GO level. Given a vector of genes, this function 
    # will return the GO profile at a specific level. 
    
    #ggo <- groupGO(gene     = genes_upregulated$ensemblid,
    #               OrgDb    = org.Mm.eg.db,
    #               ont      = "CC",
    #               level    = 3,
    #               readable = TRUE,
    #               keyType = "ENSEMBL")
    #head(ggo@result)

    
    # GO over-representation test
    # 

    for(go_ontology in list("MF", "CC", "BP")){
        # 
        cat(paste(" - Functional Enrichment - GO Over-representation test - ", go_ontology, "\n"))
        
        

        #go_ontology<- "CC"
        ego <- enrichGO(gene          = gene_list$ensemblid,
                        universe      = names(gene_list$ensemblid),
                        OrgDb         = org.Mm.eg.db,
                        ont           = go_ontology,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable = TRUE,
                        keyType = "ENSEMBL")

        # summary
        summary_filename = paste0(output_folder, title, "_", go_ontology, "_enrichment_summary",  ".txt")
        sink(summary_filename)
        cat(
            "PvalueCutoff: ", ego@pvalueCutoff, "\n",
            "pAdjustMethod: ", ego@pAdjustMethod, "\n",
            "qvalueCutoff: ", ego@qvalueCutoff, "\n",
            "Organism: ", ego@organism, "\n",
            "Ontology: GO", ego@ontology, "\n",
            "Gene: ", ego@gene, "\n",
            "Gene2symbol: ", ego@gene2Symbol, "\n",
            sep = " "
        )
        sink()  # returns output to the console
        results_filename = paste0(output_folder, title, "_", go_ontology, "_enrichment_table",  ".csv")
        write.csv(ego@result, file=results_filename)

        
        # plots
        # 
        #enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
        #barplot(ego, drop=TRUE, showCategory=30)
        #

        barplot <- barplot(ego, drop=TRUE, showCategory=30)
        filename_barplot = paste0(output_folder, title, "_", go_ontology, "_enrichment_barplot",  ".png")
        png(filename_barplot, width=800, height=800)
        print(barplot)
        dev.off()   

        
        # check if any enrichment, otherwise quit routine
        # 
        enrichedTerms <- ego@result$p.adjust <= 0.05
        numberEnrichedTerms <- as.integer(table(enrichedTerms)["TRUE"]) # count number of enriched terms
        if(numberEnrichedTerms < 2 || is.na(numberEnrichedTerms )){
            cat(paste("- No enrichment for", go_ontology, ", skipping goplot\n"))
            #next # skip and go to next for iteration
        }else{
            goplot <- goplot(ego)
            #cat("performed goplot \n")
            filename_goplot = paste0(output_folder, title, "_", go_ontology, "_enrichment_goplot",  ".png")
            png(filename_goplot, width=800, height=800)
            print(goplot)
            dev.off()   
            #cat("wrote goplot \n")
        }
        
        dotplot <- dotplot(ego)
        filename_dotplot = paste0(output_folder, title, "_", go_ontology, "_enrichment_dotplot",  ".png")
        png(filename_dotplot, width=800, height=800)
        print(dotplot)
        dev.off()

    }
    
    
    # KEGG over-representation test
    # 
    cat(" - Functional Enrichment - KEGG over-representation test\n")
    
    genes <- bitr(gene_list$ensemblid, fromType = "ENSEMBL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Mm.eg.db)
    
    kk <- enrichKEGG(gene         = genes$ENTREZID,
                     organism     = 'mmu',
                     pvalueCutoff = 0.05)
    #head(kk)
    #browseKEGG(kk, 'mmu04080') 
    
    # add in URL to Kegg pathway to dataframe
    kk_df <- kk@result
    url <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", kk_df$ID, "/", kk_df$geneID)
    # =HYPERLINK("URL") for Excel
    kk_df$url <- paste0('=HYPERLINK("', url, '")')

    # write to file
    # 
    filename = paste0(output_folder, title, "_kegg_pathways",  ".csv")
    write.csv(kk_df, file=filename)
    
    cat(" - DONE with Functional Enrichment\n")
    
}