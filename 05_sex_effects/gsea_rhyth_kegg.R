library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)


coefs=readRDS("interaction_results_readjusted.rds") 

res_wt_mf = readRDS("res_wt_mf.rds") %>% 
  pivot_wider(id_cols = c("cluster", "gene"), 
              names_from = sex, 
              values_from = padj, 
              values_fill = 1)
coefs=left_join(coefs, res_wt_mf %>% 
                  dplyr::select(cluster,gene, M,F)) %>% 
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
         region=plyr::mapvalues(cluster, regions, names(regions)),
         region=factor(region, levels=names(regions)[!duplicated(names(regions))]),
         rel_amp_dif=f_amp-m_amp,
         phi_dif=atan2(sin(f_phi-m_phi),cos(f_phi-m_phi)) %% (2 * pi),
         phi_dif_hr=12*phi_dif/pi,
         rhyth=case_when(M<0.1 & F>0.1 ~ "Male",
                         M>0.1 & F<0.1~"Female", 
                         T ~ "Both"))



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
  
gsea_res =gseKEGG(
geneList     = ranks,
organism="mmu",
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

saveRDS(all_gsea, "gsea_rhth_amp_FvM_kegg.rds")
