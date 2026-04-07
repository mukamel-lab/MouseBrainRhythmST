library(tidyverse)
library(patchwork)
setwd("~/desp1/precast/prec_c25q25g3000/deseq_geno_int")

res_sig_young=readRDS("sex_geno_young_ds_lrt/res_sig.rds")
res_sig_old=readRDS("sex_geno_old_ds_lrt/res_sig.rds")
cluster_color=readRDS("../objects/cluster_color.rds")
count_genes=function(df) {
  
  df %>% mutate(cluster=gsub("Olfactory Tubercle", "Amygdala", cluster)) %>%
    mutate(cluster=factor(cluster, levels = rev(names(cluster_color)))) %>%
    group_by(cluster) %>%
    summarise(n=n())  %>%
    ungroup() %>%
    complete(cluster,fill = list(n=0))
  
}

gene_counts=rbind(count_genes(res_sig_old) %>%
                    mutate(Age="14 months"),
                  count_genes(res_sig_young) %>% 
                    mutate(Age="7 months")
)

gene_counts=gene_counts %>% 
  dplyr::filter(cluster %in% names(cluster_color)[1:16]) %>%
  mutate(Age=relevel(factor(Age),"7 months"))


ggplot(gene_counts) + geom_col(aes(x=cluster, y=n, fill=cluster))+
  facet_wrap(~Age)+theme_bw(base_size = 16)+
  scale_fill_manual(values=cluster_color)+
  coord_flip()+
  guides(fill="none")+ylab("# DEG")+xlab("")
