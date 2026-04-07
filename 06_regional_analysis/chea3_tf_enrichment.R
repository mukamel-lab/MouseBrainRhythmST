library(rChEA3)
library(tidyverse)

interaction_results_readjusted_y =
  readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/interaction_results_readjusted_y.rds")%>%
  dplyr::filter(fdr<0.1)

ctx_drg=
  interaction_results_readjusted_y %>% 
  dplyr::filter(grepl("^Cortex", cluster), fdr<0.1) %>% .[["gene"]] %>% unique 
res = queryChEA3(ctx_drg,
                 verbose = FALSE)

res_comb=res[1:2]

res_indv=lapply(res[3:length(res)], function(df) {
  df %>% dplyr::filter(FDR<0.05)
})

sig_TFs_fdr0.05=lapply(res_indv, function(df) {
  df$TF %>% unique()
}) %>% unlist() %>% table()

res_comb=lapply(res_comb, function(df) {
  df %>% dplyr::filter(TF %in% names(sig_TFs_fdr0.05))
})

res_comb=c(res_comb, res_indv)
res_comb= lapply(res_comb, function(df) { df %>% dplyr::select(-`Query Name`)})
openxlsx::write.xlsx(res_comb, "~/desp1/supplemental_tables/TableS10 CHEA3_Cortex_DRG_DRG_APP23_vs_NTG_7_months.xlsx")
