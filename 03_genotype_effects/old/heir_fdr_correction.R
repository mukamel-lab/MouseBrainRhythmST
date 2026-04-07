library(tidyverse)


setwd("~/desp1/precast/prec_c25q25g3000")

res_wt_aw = readRDS("deseq_rhth_int2/geno_int_new/old/res_o_aw.rds") %>%
  dplyr::filter(padj<0.1)

res_wt_aw=split(res_wt_aw, res_wt_aw$cluster) 

coefs = readRDS("deseq_rhth_int2/geno_int_new/old/coefs_o.rds")

coefs=split(coefs, coefs$cluster)

coefs_readj=lapply(names(res_wt_aw), function(clst) {
  cf=coefs[[clst]] %>% 
    dplyr::filter(!is.na(padj_int)) %>%
    dplyr::filter((padj<0.001 | gene %in% res_wt_aw[[clst]]$gene ))
  cf$fdr=p.adjust(cf$pvalue_int, method = "BH")
  cf
}) %>% bind_rows()

saveRDS(coefs_readj, "deseq_rhth_int2/geno_int_new/old/interaction_results_readjusted_o.rds")

