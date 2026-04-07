library(tidyverse)


setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_age")

genes=readRDS("app_res.rds")

gene_lists=split(genes, genes$cluster)

gene_lists=lapply(gene_lists, function(x) { 
    y=x %>% na.omit() %>% arrange(desc(log2FoldChange)) 
    y=y$log2FoldChange %>% set_names(y$gene)
    y
  }) 


GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  all_res <- fgsea::fgsea(pathways = myGO,
                          stats = gene_list,
                          minSize=15, ## minimum gene set size
                          maxSize=400, ## maximum gene set size
                          nperm=10000) %>% 
    as.data.frame()
  
  fgRes= all_res %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) + 
    theme_bw()
  
  output = list("sig_res" = fgRes, "plt" = g1, "all_res"=all_res)
  return(output)
}

GO_file = "~/desp1/precast/prec_c25q25g3000/gsea/m5.go.v2023.2.Mm.symbols.gmt"

res =pbmclapply(gene_lists, function(gene_lst) {
    GSEA(gene_lst, GO_file, pval = 0.05)})

saveRDS(res, "gsea/app_gsea_res.rds",compress = F)


