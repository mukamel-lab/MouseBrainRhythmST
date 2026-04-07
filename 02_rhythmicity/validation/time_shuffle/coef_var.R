library(pbmcapply)
library(tidyverse)
library(customfuncs)
library(Matrix)

setwd("~/desp1/precast/prec_c25q25g3000/")

data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/deseq_norm_l2cpm.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data=split(meta.data, meta.data$genotype)

wt_exp=lapply(agg_exp, function(df) {
  df[, intersect(colnames(df),meta.data$WT$sample)]
})
app_exp=lapply(agg_exp, function(df) {
  df[, intersect(colnames(df),meta.data$APP23$sample)]
})

covar=function(i) {sd(i)/mean(i)}

covar_app=pbmclapply(app_exp, function(df) {
  apply(df, 1, covar)
})

covar_wt=pbmclapply(wt_exp, function(df) {
  apply(df, 1, covar)
})

for (i in names(covar_app)) {
  df=covar_app[[i]]
  df=df[!is.na(df)]
  # df=data.frame(genotype="APP23-TG",
  #               mn=mean(df), 
  #               se=sd(df)/sqrt(length(df)))
  
  df=data.frame(genotype="APP23-TG",
                covar=df)
  covar_app[[i]]=df
}


for (i in names(covar_wt)) {
  df=covar_wt[[i]]
  df=df[!is.na(df)]
  # 
  # df=data.frame(genotype="NTG",
  #               mn=mean(df), 
  #               se=sd(df)/sqrt(length(df)))
  df=data.frame(genotype="NTG",
                covar=df)
  covar_wt[[i]]=df
}

covar_app=bind_rows(covar_app, .id="cluster")
covar_wt=bind_rows(covar_wt, .id="cluster")
covar=rbind(covar_app, covar_wt)

genotype_colors=ggsci::pal_nejm()(2)
names(genotype_colors)=c("APP23-TG", "NTG")
cluster_color = readRDS("objects/cluster_color.rds")


ggplot(covar %>% 
         mutate(cluster=factor(cluster,
                               levels = rev(names(cluster_color)))),
       aes(x=cluster, y=mn, col=genotype,
                 ymin = mn - se, ymax = mn + se))+
  geom_point(position = position_dodge(width = .7))+
geom_errorbar(position=position_dodge(width = 0.7), width=0.3, linewidth=1)+
  theme_bw(base_size = 16)+coord_flip()+
scale_color_manual(values = genotype_colors, name="Genotype")+
xlab("")+ylab("Coefficients of Variation")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


ggplot(covar %>% 
         mutate(cluster=factor(cluster,
                               levels = names(cluster_color))),
       aes(x=cluster, y=covar, col=genotype))+
  geom_boxplot(position = "dodge",outlier.shape = NA)+
  theme_bw(base_size = 16)+
  scale_color_manual(values = genotype_colors, name="Genotype")+
  xlab("")+ylab("Coefficients of Variation")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 50,hjust = 1))+
  ylim(c(NA, .5))+ggtitle("CV Across Expressed Genes")
