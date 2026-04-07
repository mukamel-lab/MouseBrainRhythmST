library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/age_int_new/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)


coefs=readRDS("interaction_results_readjusted.rds") 

res_wt_yo = readRDS("res_wt_yo.rds") %>% 
  mutate(age=plyr::mapvalues(age, c("7 months", "14 months"), c("Y", "O")))%>%
  pivot_wider(id_cols = c("cluster", "gene"), 
              names_from = age, 
              values_from = padj, 
              values_fill = 1)

coefs=left_join(coefs, res_wt_yo %>% 
                  dplyr::select(cluster,gene,Y,O)) %>% 
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
         region=plyr::mapvalues(cluster, regions, names(regions)),
         region=factor(region, levels=names(regions)[!duplicated(names(regions))]),
         rel_amp_dif=o_amp-y_amp,
         phi_dif=atan2(sin(o_phi-y_phi),cos(o_phi-y_phi)) %% (2 * pi),
         phi_dif_hr=12*phi_dif/pi,
         rhyth=case_when(Y<0.1 & O>0.1 ~ "7 months",
                         Y>0.1 & O<0.1~"14 months", 
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

saveRDS(all_gsea, "gsea_rhth_amp_OvY_kegg.rds")
