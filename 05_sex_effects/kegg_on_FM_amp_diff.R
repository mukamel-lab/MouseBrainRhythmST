library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/")

coefs=readRDS("coefs_wt.rds") 
res_wt_mf = readRDS("res_wt_mf.rds") %>% 
  pivot_wider(id_cols = c("cluster", "gene"), names_from = sex, values_from = padj, values_fill = 1)
coefs=left_join(coefs, res_wt_mf %>% dplyr::select(cluster,gene, M,F))
coefs=coefs %>% dplyr::filter(M<0.1|F<0.1)
clusters <- unique(coefs$cluster)

# Initialize a list to store GSEA results per cluster
gsea_results <- vector("list", length(clusters))
names(gsea_results) <- clusters

for (cl in clusters) {
  df_cl <- coefs %>%
    filter(cluster == cl) %>%
    mutate(rank_metric = f_amp - m_amp) %>%
    filter(!is.na(rank_metric))
  
  ranks <- df_cl$rank_metric
  names(ranks) <- df_cl$gene
  ranks <- sort(ranks, decreasing = TRUE)
  
  sym2entrez <- bitr(
    names(ranks),
    fromType   = "SYMBOL",
    toType     = "ENTREZID",
    OrgDb      = org.Mm.eg.db
  ) 
  ranks=data.frame(SYMBOL=names(ranks), ranks)
  ranks=left_join(ranks,sym2entrez)
  ranks=na.omit(ranks)
  ranks=setNames(ranks$ranks, ranks$ENTREZID)
  gsea_res =gseGO(
    geneList     = ranks,
    OrgDb        = org.Mm.eg.db,
    ont          = "BP",
    keyType      = "ENTREZID",
    minGSSize    = 5,
    pvalueCutoff = 0.2
  )
  
  # 6) Store the result
  gsea_results[[cl]] <- gsea_res
}



all_gsea= lapply(names(gsea_results), function(cl) {
    gsea_results[[cl]]@result
  
  })
names(all_gsea)=names(gsea_results)
all_gsea=bind_rows(all_gsea, .id="cluster")
