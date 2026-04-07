library(tidyverse)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth/old_young_sep/ds")
cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

fls=list.dirs(recursive = F) %>% setNames(gsub("\\./","",list.dirs(recursive = F)))

gn_cts=lapply(fls, function(x) {
  readRDS(glue::glue("{x}/par_im/gene_counts.rds"))
})
gn_cts=bind_rows(gn_cts,.id="grp") %>% 
  separate(grp, c("age", "genotype", "ds"), sep = "_") %>%
  mutate(Age=plyr::mapvalues(age, c("old", "young"), c("14 months", "7 months")),
         Genotype=plyr::mapvalues(genotype, c("app", "wt"), c("APP23-TG", "NTG")),
         cluster=gsub("Olfactory Tubercle", "Amygdala", cluster),
         cluster=factor(cluster, levels=rev(names(cluster_color))),
         Age=relevel(factor(Age), "7 months"))



geno_colors=c("#0072B5FF","#BC3C29FF") %>%
  setNames(c("NTG",
             "APP23-TG"))

plt_bar=ggplot(gn_cts , aes(x=cluster, y=n, fill=Genotype))+
  geom_col(position = position_dodge())+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank())+
  theme_bw(base_size = 16)+
  facet_wrap(~Age)+coord_flip()+ylab("Number DRG")+xlab("")+
  scale_fill_manual(name="Genotype", values=geno_colors)


dta=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/spat_plots_generic/spat_data_030B.rds")

gn_cts_dif=gn_cts %>%complete(cluster=cluster) %>%pivot_wider(names_from = genotype, 
                                  values_from = n, 
                                  values_fill = 0, 
                                    
                                  id_cols = c(cluster, Age)) %>%
  mutate(dif=app-wt) %>% dplyr::rename(anot=cluster)  

dta_old=left_join(dta, gn_cts_dif %>% dplyr::filter(Age=="14 months") ,
                  by="anot") %>%
  mutate(Age="14 months", dif=replace_na(dif,0))

dta_young=left_join(dta, gn_cts_dif %>% dplyr::filter(Age=="7 months"), 
                    by="anot") %>%
  mutate(Age="7 months", dif=replace_na(dif,0))
dta_both=rbind(dta_old, dta_young)

plt_dif=ggplot(dta_both %>% mutate(Age=relevel(factor(Age), "7 months")))+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=dif), size=.5)+
  theme_void()+
  scale_color_gradient(name="", 
                       low="gray90",
                       high="darkred", 
                       limits=c(NA,NA),
                       oob=scales::squish)+
  ggtitle("")+coord_fixed()+
  theme(plot.title=element_text(size = 12, hjust = 0.5),
        legend.title=element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(.9, "lines")
  ) + facet_wrap(~Age)

plt_bar/plt_dif
