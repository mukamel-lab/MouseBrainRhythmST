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
    dplyr::select(cluster, gene, rel_amp, padj,pvalue) %>% 
    dplyr::rename(rel_amp_app=rel_amp, padj_app=padj, pvalue_app=pvalue)
  yng=left_join(yng, yng2) %>% mutate(lfc_ramp=log2(rel_amp_app/rel_amp))
  
  old=all_coefs$WT_old_harmonic
  old2=all_coefs$APP_old_harmonic %>%
    dplyr::select(cluster, gene, rel_amp,padj,pvalue) %>%
    dplyr::rename(rel_amp_app=rel_amp,padj_app=padj, pvalue_app=pvalue)
  old=left_join(old, old2)%>% mutate(lfc_ramp=log2(rel_amp_app/rel_amp))
  
   yng=split(yng, yng$gene)
   old=split(old, old$gene)
   
  
   list(yng=yng, old=old)
})

plot_rel_amp=function(gn,pcut) {
  nrclusts=rbind(all_coefs$yng[[gn]], all_coefs$yng[[gn]])
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


plot_rel_amp("Polr3gl", 1)
