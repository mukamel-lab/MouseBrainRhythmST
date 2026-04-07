library(tidyverse)
library(patchwork)
library(scico)

setwd("~/desp1/precast/prec_c25q25g3000/figures/spat_plots_generic")

dta=readRDS("spat_data_030B.rds")
dta$cluster=as.character(dta$anot)


make_coefs_mat=function(res_dir) {
  

coefs_df=readRDS(glue::glue("../../deseq_rhth/old_young_sep/{res_dir}/par_im/coefs.rds")) %>%
  bind_rows(.id="cluster")

coefs_df=left_join(coefs_df,
                         readRDS(glue::glue("../../deseq_rhth/old_young_sep/{res_dir}/par_im/res.rds")),
                         by=c("gene", "cluster"))%>%
  mutate(rel_amp=amp/log2(baseMean+1))


coefs_df$cluster=gsub("Olfactory Tubercle", "Amygdala", coefs_df$cluster)

covar_group=gsub("_harmonic", "", res_dir)
covar_group=gsub("_", " ", covar_group)
covar_group=stringi::stri_replace_all_fixed(covar_group,
                                            c("WT", "APP","young", "old"),
                                            c("NTG", "APP23-TG","7 Mo.", "14 Mo."),
                                            vectorize_all = F)

coefs_df$covar_group=covar_group


coefs_df
}

covar_groups=expand.grid(c("WT_", "APP_"), c("young","old"))
covar_groups=paste0(covar_groups$Var1, covar_groups$Var2, "_harmonic")

all_coefs=lapply(covar_groups,make_coefs_mat )
names(all_coefs)=covar_groups
all_coefs=local({
  yng=all_coefs$WT_young_harmonic
  yng2=all_coefs$APP_young_harmonic %>%
    dplyr::select(cluster, 
                  gene,
                  rel_amp, 
                  padj,
                  pvalue,
                  phi,
                  phi_hr) %>% 
    dplyr::rename(rel_amp_app=rel_amp,
                  padj_app=padj, 
                  pvalue_app=pvalue,
                  phi_app=phi,
                  phi_hr_app=phi_hr)
  
  yng=left_join(yng, yng2) %>%
    mutate(lfc_ramp=log2(rel_amp_app/rel_amp),
           phi_dif=(12/pi)*atan2(sin(phi_app-phi),cos(phi_app-phi)),
           phi_dif2=phi_hr_app-phi)

    old=all_coefs$WT_old_harmonic
    
  old2=all_coefs$APP_old_harmonic %>%
    dplyr::select(cluster,
                  gene,
                  rel_amp, 
                  padj,
                  pvalue,
                  phi,
                  phi_hr) %>% 
    dplyr::rename(rel_amp_app=rel_amp, 
                  padj_app=padj, 
                  pvalue_app=pvalue,
                  phi_app=phi,
                  phi_hr_app=phi_hr)
  
  old=left_join(old, old2)%>% 
    mutate(lfc_ramp=log2(rel_amp_app/rel_amp),
          phi_dif=(12/pi)*atan2(sin(phi_app-phi),cos(phi_app-phi)),
          phi_dif2=phi_hr_app-phi)
   
  
   list(yng=yng, old=old)
})

young_app_lfc=local({
  readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_app/all_wald/young/res.rds")%>%
    mutate(l2expr=log2(baseMean+1))%>%
    dplyr::select(cluster,gene, l2expr,log2FoldChange)%>%
  dplyr::rename(l2fc_app_base=log2FoldChange)
  
})

old_app_lfc=local({
  readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_app/all_wald/old/res.rds") %>%
    mutate(l2expr=log2(baseMean+1))%>%
    dplyr::select(cluster,gene, l2expr,log2FoldChange) %>%
    dplyr::rename(l2fc_app_base=log2FoldChange)
  
})

all_coefs$yng=left_join(all_coefs$yng,young_app_lfc,by=c("cluster","gene"))
all_coefs$old=left_join(all_coefs$old,old_app_lfc,by=c("cluster","gene"))

all_coefs$yng=split(all_coefs$yng, all_coefs$yng$gene)
all_coefs$old=split(all_coefs$old, all_coefs$old$gene)

#######################################################

plot_phase_dif=function(gn, pcut=1) {

  nrclusts=rbind(all_coefs$yng[[gn]], all_coefs$old[[gn]])
  lmt_rang=na.omit(unique(nrclusts$phi_dif))
  #lmt_rang=max(abs(quantile(lmt_rang, c(0.01, .99))))
  
  lmt_rang=max(abs(lmt_rang), na.rm = T)
  
  nrclusts= nrclusts %>%
    dplyr::filter(padj>pcut, padj_app>pcut)
  
  nrclusts=unique(nrclusts$cluster)
  
  
  df=left_join(dta,all_coefs$yng[[gn]] )%>%
   mutate(phi_dif=replace_na(phi_dif, 0)) %>%
    mutate(phi_dif=ifelse(cluster %in% nrclusts, NA, phi_dif))


p_young=ggplot(df)+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=phi_dif), size=0.5)+
  theme_void()+
  scale_color_scico(palette = 'vik', name="Phase\nDifference", midpoint = 0,
                    limits=c(-1*lmt_rang, lmt_rang),
                    oob=scales::squish) +
  ggtitle("7 Mo.")+coord_fixed()+
  theme(plot.title=element_text(size = 12, hjust = 0.5),
        legend.title=element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(.9, "lines")
        )


df=left_join(dta,all_coefs$old[[gn]] ) %>%
  mutate(phi_dif=replace_na(phi_dif, 0))%>%
  mutate(phi_dif=ifelse(cluster %in% nrclusts, NA, phi_dif))
p_old=ggplot(df)+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=phi_dif), size=0.5)+
  theme_void()+
  scale_color_scico(palette = 'vik', name="Phase\nDifference", midpoint = 0,
                    limits=c(-1*lmt_rang, lmt_rang),
                    oob=scales::squish) +
  ggtitle("14 Mo.")+coord_fixed()+
  theme(plot.title=element_text(size = 12, hjust = 0.5),
        legend.title=element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(.9, "lines")
  )

(p_young | p_old) +plot_layout(guides = "collect")+
  plot_annotation(title = gn,
                  theme = theme(plot.title =
                                  element_text(size = 16, 
                                               hjust = 0.5)))
}

################################################
plot_rel_amp=function(gn,pcut) {
  nrclusts=rbind(all_coefs$yng[[gn]], all_coefs$old[[gn]])
  lmt_rang=na.omit(unique(nrclusts$lfc_ramp))
  #lmt_rang=max(abs(quantile(lmt_rang, c(0.1, .9))))
  
  lmt_rang=max(abs(lmt_rang))
  
  nrclusts= nrclusts %>%
    dplyr::filter(padj>pcut, padj_app>pcut)
  
  
  nrclusts=unique(nrclusts$cluster)
  df=left_join(dta,all_coefs$yng[[gn]] )%>%
    mutate(lfc_ramp=replace_na(lfc_ramp, 0)) %>%
    mutate(lfc_ramp=ifelse(cluster %in% nrclusts, NA, lfc_ramp))
  
  
  p_young=ggplot(df)+
    geom_point(aes(y=-1*imagerow,x=imagecol, col=lfc_ramp), size=0.5)+
    theme_void()+
    scale_color_scico(palette = 'bam', name="Log2FC\nRel. Amp.", midpoint = 0,
                      limits=c(-1*lmt_rang, lmt_rang),
                      oob=scales::squish) +
    ggtitle("7 Mo.")+coord_fixed()+
    theme(plot.title=element_text(size = 12, hjust = 0.5),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(.9, "lines")
    )
  
  
  df=left_join(dta,all_coefs$old[[gn]] ) %>%
    mutate(lfc_ramp=replace_na(lfc_ramp, 0))%>%
    mutate(lfc_ramp=ifelse(cluster %in% nrclusts, NA, lfc_ramp))
  
  
  p_old=ggplot(df)+
    geom_point(aes(y=-1*imagerow,x=imagecol, col=lfc_ramp), size=0.5)+
    theme_void()+
    scale_color_scico(palette = 'bam', name="Log2FC\nRel. Amp.", midpoint = 0,
                      limits=c(-1*lmt_rang, lmt_rang),
                      oob=scales::squish) +
    ggtitle("14 Mo.")+coord_fixed()+
    theme(plot.title=element_text(size = 12, hjust = 0.5),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(.9, "lines")
    )
  
  (p_young | p_old) +plot_layout(guides = "collect")+
    plot_annotation(title = gn,
                    theme = theme(plot.title =
                                    element_text(size = 16, 
                                                 hjust = 0.5)))
}

##########################################################

plot_lfc=function(gn,pcut) {
  nrclusts=rbind(all_coefs$yng[[gn]], all_coefs$old[[gn]])
  lmt_rang=na.omit(unique(nrclusts$l2fc_app_base))
  #lmt_rang=max(abs(quantile(lmt_rang, c(0.1, .9))))
  max_l2expr=max(nrclusts$l2expr, na.rm = T)
  lmt_rang=max(abs(lmt_rang))
  
  nrclusts= nrclusts %>%
    dplyr::filter(padj>pcut, padj_app>pcut)
  
  
  nrclusts=unique(nrclusts$cluster)
  df=left_join(dta,all_coefs$yng[[gn]] )%>%
    mutate(lfc_ramp=replace_na(lfc_ramp, 0)) %>%
    mutate(lfc_ramp=ifelse(cluster %in% nrclusts, NA, lfc_ramp))
  
  
  p_young=ggplot(df)+
    geom_point(aes(y=-1*imagerow,x=imagecol, col=l2fc_app_base, size=l2expr))+
    theme_void()+
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name="Log2FC\nAPP23 vs. NTG", 
                         limits=c(-1*lmt_rang, lmt_rang),
                         oob=scales::squish) +
    ggtitle("7 Mo.")+coord_fixed()+
    theme(plot.title=element_text(size = 12, hjust = 0.5),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(.9, "lines")
    ) +
    scale_size(range = c(0,.5),
               name = "log2(Mean CPM)",
               limits = c(0,max_l2expr+0.1)
                )

  
  
  df=left_join(dta,all_coefs$old[[gn]] ) %>%
    mutate(lfc_ramp=replace_na(lfc_ramp, 0))%>%
    mutate(lfc_ramp=ifelse(cluster %in% nrclusts, NA, lfc_ramp))
  
  
  p_old=ggplot(df)+
    geom_point(aes(y=-1*imagerow,x=imagecol, col=l2fc_app_base, size=l2expr))+
    theme_void()+
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         name="Log2FC\nAPP23 vs. NTG", 
                         limits=c(-1*lmt_rang, lmt_rang),
                         oob=scales::squish) +
    
    ggtitle("14 Mo.")+coord_fixed()+
    theme(plot.title=element_text(size = 12, hjust = 0.5),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(.9, "lines")
    )+
            scale_size(range = c(0,0.5),
                       name = "log2(Mean CPM)",
                       limits = c(0,max_l2expr+0.1))
  
  (p_young | p_old) +plot_layout(guides = "collect")+
    plot_annotation(title = gn,
                    theme = theme(plot.title =
                                    element_text(size = 16, 
                                                 hjust = 0.5), legend.box = "horizontal"))
}

(plot_phase_dif("Atf3", 1)/
plot_rel_amp("Atf3", 1)) +plot_annotation(title = "Atf3",
                                         theme = theme(plot.title =
                                                         element_text(size = 16, 
                                                                      hjust = 0.5)))

gn="Per2"
((plot_rel_amp(gn, 1)/
(plot_phase_dif(gn, 1)&ggtitle(""))/

(plot_lfc(gn, 1)&ggtitle(""))
)&theme(legend.position = "left"))+
  plot_annotation(title = gn,
                  theme = theme(plot.title =
                                  element_text(size = 16, 
                                               hjust = 0.5)))


plot_rel_amp("Foxo3", 1)
plot_phase_dif("Midn", 1)


plot_rel_amp("Foxo3", 1)/
  plot_rel_amp("Mical3", 1)

