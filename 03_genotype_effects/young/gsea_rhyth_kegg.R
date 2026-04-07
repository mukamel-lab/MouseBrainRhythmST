library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)

res_wt_aw = readRDS("old/res_o_aw.rds") %>%
  dplyr::filter(padj<0.05)

res_wt_aw=split(res_wt_aw, res_wt_aw$cluster) 

coefs = readRDS("old/coefs_o.rds") 

coefs=split(coefs, coefs$cluster)

coefs=lapply(names(res_wt_aw), function(clst) {
  cf=coefs[[clst]] %>% 
    dplyr::filter((padj<0.05|gene %in% res_wt_aw[[clst]]$gene ))
  cf
}) %>% bind_rows()

coefs=coefs %>%
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
         region=plyr::mapvalues(cluster, regions, names(regions)),
         region=factor(region, levels=names(regions)[!duplicated(names(regions))]),
         rel_amp_dif=app_amp-wt_amp)



clusters = unique(coefs$cluster)

gsea_results = vector("list", length(clusters))
names(gsea_results) = clusters

for (cl in clusters) {
  set.seed(1)
  df_cl = coefs %>%
    filter(cluster == cl) %>%
    mutate(rank_metric = rel_amp_dif) %>%
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
  
# gsea_res =gseKEGG(
# geneList     = ranks,
# organism="mmu",
# minGSSize    = 10,
# eps=0,
# pvalueCutoff = 1,
# seed=T,
# pAdjustMethod = "BH",
# 
# )

gsea_res =gseGO(
  geneList     = ranks,
  OrgDb        = org.Mm.eg.db,
  ont          = "ALL",
  keyType      = "ENTREZID",
  minGSSize    = 10,
  pvalueCutoff = 0.1,
  eps = 0,
  seed=T,
  pAdjustMethod = "BH",

)

  
  
  gsea_results[[cl]] = gsea_res
}


all_gsea = lapply(gsea_results, function(cl) {
  cl@result
}) %>% bind_rows(.id = "cluster")

saveRDS(all_gsea, "old/old_gsea_rhth_amp_app_vs_wt_GO.rds")

