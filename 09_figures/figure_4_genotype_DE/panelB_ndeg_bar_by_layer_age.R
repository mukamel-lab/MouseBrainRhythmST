library(tidyverse)
library(patchwork)

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_app/ds_bootstrap")

deg_res=list(young=lapply(list.files("7_months_objects/", full.names = T), readRDS) %>%
  setNames(gsub("\\.rds", "",list.files("7_months_objects/"))) %>% 
             bind_rows(.id="run" ),
old=lapply(list.files("14_months_objects/", full.names = T), readRDS) %>%
  setNames(gsub("\\.rds", "",list.files("14_months_objects/"))) %>% 
  bind_rows(.id="run" )) %>% bind_rows(.id="Age") %>%
  mutate(Age=plyr::mapvalues(Age, c("young", "old"), paste0(c(7, 14), " months")))
deg_res=deg_res %>% group_by(Age, cluster, run) %>% tally() %>%
  ungroup() %>% 
  group_by(Age, cluster) %>%
  summarise(mean=mean(n), se=sd(n)/sqrt(n()))


deg_res$cluster=gsub(" \\(L6b, CLA, EP\\)", "",
                 gsub("Dentate Gyrus",  "DG",
                      gsub("Cortex Layer ", 
                           "CTX L", deg_res$cluster)))

age_colors=c(RColorBrewer::brewer.pal(12, "Paired")[11],
             RColorBrewer::brewer.pal(3, "BrBG")[1])

ggplot(deg_res %>% mutate(), aes(x=cluster, y=mean, fill=Age,
                    ymin = mean - se, ymax = mean + se))+
  geom_col(position=position_dodge(width = 0.6), width = 0.6)+
  geom_errorbar(position=position_dodge(width = 0.6), width=0.3 ) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank())+
  theme_bw(base_size = 16)+theme(axis.text.x = element_blank(),
                                 axis.ticks.x = element_blank())+
  scale_fill_manual(name="Age", values=age_colors)+xlab("")+ylab("")





