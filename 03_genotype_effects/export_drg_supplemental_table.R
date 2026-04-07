library(tidyverse)

interaction_results_readjusted_y <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds")%>% dplyr::filter(fdr<0.1)%>%
  dplyr::select(-padj_int)
interaction_results_readjusted_o <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/old/interaction_results_readjusted_o.rds")%>% dplyr::filter(fdr<0.1)%>%
  dplyr::select(-padj_int)

openxlsx::write.xlsx(list("7 months"= interaction_results_readjusted_y, "14 months"=interaction_results_readjusted_o), 
                     "~/desp1/supplemental_tables/APP23_vs_NTG_DRG.xlsx")
