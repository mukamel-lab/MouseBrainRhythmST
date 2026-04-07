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

res_both_sig=res_both %>% group_by(cluster) %>% dplyr::filter((app_padj<0.05| wt_padj<0.05))
res_both_sig$sig="a"
res_both_sig$sig[res_both_sig$app_padj<0.05]="APP23-TG" 
res_both_sig$sig[res_both_sig$wt_padj<0.05]="NTG"
res_both_sig$sig[res_both_sig$app_padj<0.05 & res_both_sig$wt_padj<0.05]="Both" 

res_both_sig$sig=factor(res_both_sig$sig, levels=unique(res_both_sig$sig)[c(1,3,2)])

res_both_sig$relamp_lfc=log2(res_both_sig$app_rel_amp/res_both_sig$wt_rel_amp)
res_both_sig  %>% dplyr::filter(sig=="Both") %>%
  mutate(cluster=factor(cluster, levels=rev(cluster_order))) %>% ggplot(aes(x=cluster, y=relamp_lfc, fill=cluster))+
  scale_fill_manual(values=cluster_color)+
  geom_violin()+geom_boxplot(width=0.2) +
  coord_flip()+theme_bw()+ylim(c(-1,1))+guides(fill="none")+
  theme(axis.text = element_text(size = 16))
