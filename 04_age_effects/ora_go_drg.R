library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)

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
  
  ##########################################################
  df=coefs %>%
    filter(cluster == cl)
  
  df$gene2=limma::alias2SymbolTable(df$gene, species = "Mm")
  
  df=df %>% dplyr::filter(!is.na(gene2))
  
  drg= df %>%
    filter(fdr<0.1) %>% pull(gene2)
  
  univ=df$gene2
  
  
  sym2entrez = bitr(
    df$gene2,
    fromType   = "SYMBOL",
    toType     = "ENTREZID",
    OrgDb      = org.Mm.eg.db
  )
  
  drg=plyr::mapvalues(drg, sym2entrez$SYMBOL, sym2entrez$ENTREZID)
  univ=plyr::mapvalues(univ, sym2entrez$SYMBOL, sym2entrez$ENTREZID)
  
  
  
  ego = enrichGO(
    gene          = drg,
    universe      = univ,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "ALL",
    pAdjustMethod = "BH",
    pool=T
  )
  
  
  
  
  ekegg = enrichKEGG(
    gene          = drg,
    universe      = univ, 
    organism = "mmu"
  )
  
  ereact= enrichPathway(gene=drg, 
                        pvalueCutoff = 0.05,
                        readable=TRUE, 
                        organism = "mouse",
                        universe = univ
  )
  
  
  gsea_results[[cl]] =list(go=ego, kegg=ekegg, react=ereact)
}

saveRDS(gsea_results,"ora_drg.rds")

gsea_comb=lapply(gsea_results, function(x) {
  lapply(x, function(x) {
    if(is.null(x)) {return(NULL)}
    x@result
  }) %>% bind_rows(.id="pthw")
}) %>% bind_rows(.id="cluster")

saveRDS(gsea_comb,"ora_drg_dfs.rds")
