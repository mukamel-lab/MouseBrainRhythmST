library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)  
library(openxlsx)
library(tidyverse)
OrgDb = org.Mm.eg.db

setwd("~/desp1/precast/prec_c25q25g3000/deseq_reg_int2/")

res = readRDS("regional_rhythms_APP23.rds")

res_kegg = lapply(res, function(df) {
  z = df
  
  uni_map = bitr(unique(z$gene),
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = OrgDb)
  
  z$eid=plyr::mapvalues(z$gene, uni_map$SYMBOL,uni_map$ENTREZID)
  z=z %>% arrange(desc(amp))

  geneList=z$amp
  names(geneList)=z$eid
  
  gk = gseKEGG(
    geneList     = geneList,
    organism     = "mmu",
    pvalueCutoff = 0.2,
    scoreType = "pos"
  )
  
  gk@result %>% as.data.frame() 
}) %>% bind_rows(.id = "cluster")

## GO GSEA
res_go = lapply(res, function(df) {
  z = df
  
  uni_map = bitr(unique(z$gene),
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = OrgDb)
  
  z$eid=plyr::mapvalues(z$gene, uni_map$SYMBOL,uni_map$ENTREZID)
  z=z %>% arrange(desc(amp))
  
  geneList=z$amp
  names(geneList)=z$eid
  
  
  ggo_BP = gseGO(
    geneList     = geneList,
    OrgDb        = OrgDb,
    keyType      = "ENTREZID",
    ont          = "BP",
    pvalueCutoff = 0.2,
    scoreType = "pos"
  )
  
  ggo_CC = gseGO(
    geneList     = geneList,
    OrgDb        = OrgDb,
    keyType      = "ENTREZID",
    ont          = "CC",
    pvalueCutoff = 0.2,
    scoreType = "pos"
  )
  
  ggo_MF = gseGO(
    geneList     = geneList,
    OrgDb        = OrgDb,
    keyType      = "ENTREZID",
    ont          = "MF",
    pvalueCutoff = 0.2,
    scoreType = "pos"
  )
  
  ego = lapply(
    list("MF" = ggo_MF, "CC" = ggo_CC, "BP" = ggo_BP),
    function(df) df@result %>% as.data.frame()
  ) %>% bind_rows(.id = "ontology")
  
}) %>% bind_rows(.id = "cluster")

saveRDS(res_kegg, "gsea_res/gsea_res_kegg_APP23.rds")
saveRDS(res_go,   "gsea_res/gsea_res_go_APP23.rds")
