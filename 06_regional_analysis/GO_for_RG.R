library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)  
library(openxlsx)
library(tidyverse)
OrgDb = org.Mm.eg.db

setwd("~/desp1/precast/prec_c25q25g3000/deseq_reg_int2/")

res = readRDS("regional_rhythms_NTG.rds")


res_kegg=lapply(res, function(df) {
  z=df
  
  z1=df %>% dplyr::filter(padj<0.1)
  
  
  gene_symbols = unique(z1$gene)      
  all_symbols  = unique(z$gene)     
  
  fg_map  = bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
  uni_map = bitr(all_symbols,  fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
  
  fg_entrez  = unique(fg_map$ENTREZID)
  uni_entrez = unique(uni_map$ENTREZID)
  
  
  ek = enrichKEGG(
    gene         = fg_entrez,
    organism     = "mmu",       # "mmu" for mouse, "hsa" for human
    universe     = uni_entrez,
    pvalueCutoff = 0.2,
    qvalueCutoff = .2
  )
  
  
  ek@result%>% as.data.frame() 
}) %>% bind_rows(.id="cluster")

res_go=lapply(res, function(df) {
  z=df
  
  z1=df %>% dplyr::filter(padj<0.1)
  
  
  gene_symbols = unique(z1$gene)      
  all_symbols  = unique(z$gene)     
  
  fg_map  = bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
  uni_map = bitr(all_symbols,  fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
  
  fg_entrez  = unique(fg_map$ENTREZID)
  uni_entrez = unique(uni_map$ENTREZID)
  
  
  ego_BP <- enrichGO(
    gene          = fg_entrez,
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = "BP",
    universe      = uni_entrez,
    pvalueCutoff  = 0.2,
    qvalueCutoff  = 0.2
  )
  
  
  ego_CC <- enrichGO(
    gene          = fg_entrez,
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = "CC",
    universe      = uni_entrez,
    pvalueCutoff  = 0.2,
    qvalueCutoff  = 0.2
  )
  
  
  ego_MF <- enrichGO(
    gene          = fg_entrez,
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = "MF",
    universe      = uni_entrez,
    pvalueCutoff  = 0.2,
    qvalueCutoff  = 0.2
  )
  ego=lapply(list("MF"= ego_MF, "CC"=ego_CC, "BP"=ego_BP), function(df){
    df@result%>% as.data.frame() 
  }
    ) %>% bind_rows(.id="ontology")
  
}) %>% bind_rows(.id="cluster")

saveRDS(res_kegg, "enrich_res_kegg_NTG.rds")
saveRDS(res_go, "enrich_res_go_NTG.rds")