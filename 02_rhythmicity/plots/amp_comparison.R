library(tidyverse)
library(patchwork)
setwd("/home/agelber/desp1/precast/prec_c17q17g3000/deseq")
coefs_wt=readRDS("coefs_wt.rds")
coefs_app=readRDS("coefs_app.rds")

coefs_all=lapply(names(coefs_wt), function(cl){
  x=coefs_wt[[cl]]%>% mutate(rel_amp=amp/(Intercept))
  y=coefs_app[[cl]] %>% mutate(rel_amp=amp/(Intercept))
  left_join(x,y, by="gene", suffix=c(".wt", ".app"))
})

names(coefs_all)=names(coefs_wt)
sig_both_g=readRDS("sig_g_both.rds")
plts=lapply(names(coefs_all), function(cl){
  
  ggplot(coefs_all[[cl]] %>% dplyr::filter(gene %in% sig_both_g)) +
  geom_point(aes(x=rel_amp.wt, y=rel_amp.app),size=0.4)+
  xlim(0,1)+ylim(0,1)+coord_fixed()+theme_bw()+xlab("Ntg Relative Amplitude")+
    ylab("TG-APP23 Relative Amplitude")+ggtitle(cl)+
    geom_abline(intercept = 0, slope = 1)
})
names(plts)=names(coefs_all)
(plts$`Dentate Gyrus` | plts$`Cortex Layer 5/6a`)/
  (plts$`Cortex Layer 2/3` | plts$`Piriform Cortex`)
