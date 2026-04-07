library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)

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



  set.seed(1)
  
  
  ##########################################################
  df=coefs %>%
    filter(cluster == "Piriform")
  
  df$gene2=limma::alias2SymbolTable(df$gene, species = "Mm")
  
  df=df %>% dplyr::filter(!is.na(gene2))
  
  drg= df %>%
    filter(fdr<0.1, rhyth=="Male") %>% pull(gene2)
  
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


saveRDS(list(go=ego, kegg=ekegg, react=ereact), "ora_male_drg_piriform.rds")


#############################################################################


df=coefs %>%
  filter(cluster == "Piriform")

df$gene2=limma::alias2SymbolTable(df$gene, species = "Mm")

df=df %>% dplyr::filter(!is.na(gene2))

drg= df %>%
  filter(fdr<0.1, rhyth!="Male") %>% pull(gene2)

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


saveRDS(list(go=ego, kegg=ekegg, react=ereact), "ora_female_drg_piriform.rds")
