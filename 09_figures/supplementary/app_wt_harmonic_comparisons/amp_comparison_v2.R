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

res_both_sig_split=split(res_both_sig, res_both_sig$cluster)


#########
#plts
###########

plts=lapply(names(res_both_sig_split), function(cl) {
  df=res_both_sig_split[[cl]]
  ggplot(df)+geom_point(aes(x=wt_rel_amp,y=app_rel_amp, col=sig), size=.2)+
    coord_fixed()+
    theme_bw()+
    geom_abline(intercept = 0, slope=1)+
    ggtitle(cl)+
    xlab("")+
    ylab("")+
    ggsci::scale_color_npg(name="FDR<0.05")+
    theme(panel.grid.minor = element_blank())+
    xlim(c(NA, max(c(df$wt_rel_amp,df$app_rel_amp))))+
    ylim(c(NA, max(c(df$wt_rel_amp,df$app_rel_amp))))+scale_x_log10()+scale_y_log10()
  
  
})

names(plts)=names(res_both_sig_split)
names(plts)[names(plts)=="Olfactory Tubercle"]="Amgydala"
plts=plts[cluster_order]
wrap_plots(plts[1:6], ncol = 2)+plot_layout(guides="collect")
