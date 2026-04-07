library(tidyverse)
library(patchwork)

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


coefs_df=split(coefs_df,coefs_df$gene)
coefs_df
}

covar_groups=expand.grid(c("WT_", "APP_"), c("young","old"))
covar_groups=paste0(covar_groups$Var1, covar_groups$Var2, "_harmonic")

all_coefs=lapply(covar_groups,make_coefs_mat )

plot_rel_amp=function(gn, coef_df) {
  
  df=left_join(dta,coef_df[[gn]] ) %>%
    mutate(rel_amp=replace_na(rel_amp, 0))

ggplot(df)+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=log(rel_amp+1)), size=0.5)+
  theme_void()+
  scale_color_gradient(name="", 
                       low="gray90",
                       high="darkred", 
                       limits=c(0,quantile(unique(df$rel_amp), .95)),
                       oob=scales::squish)+
  ggtitle(df$covar_group[1])+coord_fixed()+
  theme(plot.title=element_text(size = 12, hjust = 0.5),
        legend.title=element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(.9, "lines")
        )

}

gn="Cd4"
plts=lapply(all_coefs,function(cf_df) {plot_rel_amp(gn, cf_df)})

wrap_plots(plts, nrow = 2)+
  plot_annotation(title = glue::glue("{gn} Relative Amplitude"),
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
)
