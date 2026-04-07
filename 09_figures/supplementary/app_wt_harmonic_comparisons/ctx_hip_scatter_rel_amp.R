library(tidyverse)

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/")


cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

res_sig_app = left_join(
  readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/res.rds"),
  bind_rows( readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/coefs.rds"), .id="cluster"),
  by=c("gene", "cluster")) %>% mutate(df="app",rel_amp=amp/log2(baseMean+1))

res_sig_wt=left_join(
  readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res.rds"),
  bind_rows( readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/coefs.rds"), .id="cluster"),
  by=c("gene", "cluster"))%>% mutate(df="wt", rel_amp=amp/log2(baseMean+1))

gene_df=bind_rows(res_sig_app, res_sig_wt)

gene_df=gene_df %>% pivot_wider(id_cols = c(gene, cluster),
                                names_from = df, 
                                values_from = colnames(gene_df)[c(3, 8, 15, 16,18)]) %>%
  mutate(sig=case_when(padj_app<0.01 & padj_wt<0.01 ~ "Both",
                       padj_app<0.01 & padj_wt>0.01 ~ "APP23-TG",
                       padj_app>0.01 & padj_wt<0.01 ~ "NTG",
                       padj_app>0.01 & padj_wt>0.01 ~ "None")) %>%
  dplyr::filter(sig!="None", (rel_amp_wt>0.05 | rel_amp_app>0.05))

gene_df=split(gene_df, gene_df$cluster)

ctx=bind_rows(gene_df[names(cluster_color)[c(1:11)]], .id = "cluster")

ctx=ctx %>% dplyr::filter(sig %in% c("Both"))

ggplot(ctx %>% 
         mutate(cluster=factor(cluster, levels =names(cluster_color)[c(1:11)] )), 
       aes(x = rel_amp_wt, y = rel_amp_app,col = cluster)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = cluster_color,name="")+
   theme_bw(base_size = 14)+
  coord_fixed()+
  scale_x_log10() + scale_y_log10()+
  xlab("log10(Rel. Amp. NTG)")+
  ylab("log10(Rel. Amp. APP23-TG)")+
  ggtitle("Relative Amplitude for Shared RG")+
  guides(col= guide_legend(nrow = 4)) +
  theme(legend.position = "bottom")
