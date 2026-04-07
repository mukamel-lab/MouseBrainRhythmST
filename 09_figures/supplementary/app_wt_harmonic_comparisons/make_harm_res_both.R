library(tidyverse)
library(patchwork)
######
# APP
######
res <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/res.rds")
coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/coefs.rds") %>%
  bind_rows(.id="cluster")


res_app=left_join(res, coefs %>% 
                    dplyr::select(cluster,gene, phi,phi_hr,amp)) %>%
  mutate(rel_amp=amp/log2(baseMean+1))

colnames(res_app)[3:ncol(res_app)]=paste0("app_", colnames(res_app)[3:ncol(res_app)])


######
# WT
######
res <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res.rds")
coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/coefs.rds") %>%
  bind_rows(.id="cluster")


res_WT=left_join(res, coefs %>% dplyr::select(cluster,gene, phi,phi_hr,amp)) %>%
  mutate(rel_amp=amp/log2(baseMean+1))

colnames(res_WT)[3:ncol(res_WT)]=paste0("wt_", colnames(res_WT)[3:ncol(res_WT)])

#####
#Both
######

res_both=left_join(res_app, res_WT)

sig_level=0.001
res_both_sig=res_both %>% group_by(cluster) %>% dplyr::filter((app_padj<sig_level| wt_padj<sig_level))
res_both_sig$sig="a"
res_both_sig$sig[res_both_sig$app_padj<sig_level]="APP23-TG" 
res_both_sig$sig[res_both_sig$wt_padj<sig_level]="NTG"
res_both_sig$sig[res_both_sig$app_padj<sig_level & res_both_sig$wt_padj<sig_level]="Both" 

res_both_sig$sig=factor(res_both_sig$sig, levels=unique(res_both_sig$sig)[c(1,3,2)])

res_both_sig_split=split(res_both_sig, res_both_sig$cluster)
saveRDS(res_both_sig_split$`Dentate Gyrus-sg`, "~/desp1/precast/prec_c25q25g3000/figures/fig3_app_wt_harm/both_sig_dgsg.rds")
