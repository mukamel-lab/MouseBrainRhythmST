library(tidyverse)
setwd("~/desp1/precast/prec_c25q25g3000")

res_wt_aw = readRDS("deseq_rhth_int2/geno_int_new/young/res_y_aw.rds") %>%
  dplyr::filter(padj<0.1) %>% group_by(gene, cluster) %>% dplyr::filter(n()==2)

res_wt_aw=split(res_wt_aw, res_wt_aw$cluster) 

coefs = readRDS("deseq_rhth_int2/geno_int_new/young/coefs_y.rds")

coefs=split(coefs, coefs$cluster)

coefs_readj=lapply(names(res_wt_aw), function(clst) {
  res_wt_aw[[clst]]
  cf=coefs[[clst]] %>% 
    dplyr::filter(gene %in% res_wt_aw[[clst]]$gene ) %>% ungroup()
  cf$fdr=p.adjust(cf$pvalue_int, method = "BH")
  cf
}) %>% bind_rows()

saveRDS(coefs_readj, "deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y_doub_rhyth.rds")

