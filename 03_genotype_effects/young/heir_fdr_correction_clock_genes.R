library(tidyverse)
setwd("~/desp1/precast/prec_c25q25g3000")

clock_g=readRDS("~/desp1/precast/precast_final_with_ros_caud/kegg_clock_genes_all.rds")


coefs = readRDS("deseq_rhth_int2/geno_int_new/young/coefs_y.rds")

coefs=split(coefs, coefs$cluster)

coefs_readj=lapply(names(res_wt_aw), function(clst) {
  cf=coefs[[clst]] %>% 
    dplyr::filter(gene %in% clock_g)
  cf$fdr=p.adjust(cf$pvalue_int, method = "BH")
  cf
}) %>% bind_rows()

saveRDS(coefs_readj, "deseq_rhth_int2/geno_int_new/young/interaction_results_clock_g_y.rds")

