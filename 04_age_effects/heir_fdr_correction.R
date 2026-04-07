library(tidyverse)
library(ggh4x)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000")

res_wt_yo = readRDS("deseq_rhth_int2/age_int_new/res_wt_yo.rds") %>%
  dplyr::filter(padj<0.1)

res_wt_yo=split(res_wt_yo, res_wt_yo$cluster) 

coefs = readRDS("deseq_rhth_int2/age_int_new/coefs_wt.rds")

coefs=split(coefs, coefs$cluster)

coefs_readj=lapply(names(res_wt_yo), function(clst) {
  cf=coefs[[clst]] %>% 
    dplyr::filter(( gene %in% res_wt_yo[[clst]]$gene ))
  cf$fdr=p.adjust(cf$pvalue_int, method = "BH")
  cf
}) %>% bind_rows()

saveRDS(coefs_readj, "deseq_rhth_int2/age_int_new/interaction_results_readjusted.rds")


