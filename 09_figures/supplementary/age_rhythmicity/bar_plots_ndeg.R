library(tidyverse)
library(patchwork)
setwd("~/desp1/precast/prec_c25q25g3000/deseq_geno_age_int/")

res_sig=readRDS("age_geno_ds25_lrt/res_sig.rds")

cluster_color=readRDS("../objects/cluster_color.rds")
count_genes=function(df) {
  
  df %>% mutate(cluster=gsub("Olfactory Tubercle", "Amygdala", cluster)) %>%
    mutate(cluster=factor(cluster, levels = rev(names(cluster_color)))) %>%
    group_by(cluster) %>%
    summarise(n=n())  %>%
    ungroup() %>%
    complete(cluster,fill = list(n=0))
  
}

gene_counts=count_genes(res_sig) 

gene_counts=gene_counts %>% 
  dplyr::filter(cluster %in% names(cluster_color)[1:16]) %>%
  mutate(Age=relevel(factor(Age),"7 months"))


ggplot(gene_counts) + geom_col(aes(x=cluster, y=n, fill=cluster))+
  facet_wrap(~Age)+theme_bw(base_size = 16)+
  scale_fill_manual(values=cluster_color)+
  coord_flip()+
  guides(fill="none")+ylab("# DEG")+xlab("")
