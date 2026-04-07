library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/age_int_new/main_effect/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color) = plyr::mapvalues(names(cluster_color), short_names, names(short_names))

degs_full = readRDS("coefs_wt_main_effect_age.rds") %>%
  mutate(cluster = plyr::mapvalues(cluster, short_names, names(short_names)))

clusters = unique(degs_full$cluster)

gsea_results = vector("list", length(clusters))
names(gsea_results) = clusters
gsea_results_simplified=gsea_results

for (cl in clusters) {
  set.seed(1)
  df_cl = degs_full %>%
    filter(cluster == cl) %>%
    mutate(rank_metric = log2FoldChange) %>%
    filter(!is.na(rank_metric))
  
  
  ranks = df_cl$rank_metric
  
  names(ranks) = df_cl$gene
  ranks = sort(ranks, decreasing = TRUE)
  
  sym2entrez <- bitr(
    names(ranks),
    fromType   = "SYMBOL",
    toType     = "ENTREZID",
    OrgDb      = org.Mm.eg.db
  )
  ranks = data.frame(SYMBOL = names(ranks), ranks)
  ranks = left_join(ranks, sym2entrez)
  ranks = na.omit(ranks)
  ranks = setNames(ranks$ranks, ranks$ENTREZID)
  
gsea_res = gseKEGG(
geneList     = ranks,
organism     = "mmu",
minGSSize    = 10,
pvalueCutoff = 0.1,
verbose      = FALSE,
seed = T,
eps=0
)



  
  gsea_results[[cl]] = gsea_res

}


all_gsea = lapply(gsea_results, function(cl) {
  cl@result
}) %>% bind_rows(.id = "cluster")

saveRDS(all_gsea, "gsea_main_effect_age_kegg.rds")
