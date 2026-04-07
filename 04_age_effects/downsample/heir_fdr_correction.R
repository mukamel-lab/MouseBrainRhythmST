library(tidyverse)
library(ggh4x)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000")

res_wt_mf = readRDS("deseq_rhth_int2/sex_int_geno_spec/downsample/res_wt_mf_ds.rds") %>%
  dplyr::filter(padj<0.1)

res_wt_mf=split(res_wt_mf, res_wt_mf$cluster) 

coefs = readRDS("deseq_rhth_int2/sex_int_geno_spec/downsample/coefs_wt_ds.rds")

coefs=split(coefs, coefs$cluster)

coefs_readj=lapply(names(res_wt_mf), function(clst) {
  cf=coefs[[clst]] %>% 
    dplyr::filter(!is.na(padj_int)) %>%
    dplyr::filter((padj<0.1 | gene %in% res_wt_mf[[clst]]$gene ))
  cf$fdr=p.adjust(cf$pvalue_int, method = "BH")
  cf
}) %>% bind_rows()

saveRDS(coefs_readj, "deseq_rhth_int2/sex_int_geno_spec/downsample/interaction_results_readjusted_ds.rds")


