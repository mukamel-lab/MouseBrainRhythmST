library(tidyverse)

cluster_color <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2")
sig_genes=vroom::vroom("../deseq_rhth_int/for_share/full/results.csv.gz") %>%
  filter(grepl("Cortex Layer 2/3", cluster), padj<0.1) 

sig_genes=sig_genes$gene

wt=rbind(readRDS("WT_young_harmonic/par_im/coefs.rds") %>%
  bind_rows(.id="cluster") %>%
  filter(grepl("Cortex L", cluster), gene %in% sig_genes) %>%
  select(cluster,gene,phi_hr) %>% 
  mutate(Age="7 months") %>% dplyr::rename(wt=phi_hr),

readRDS("WT_old_harmonic/par_im/coefs.rds") %>%
  bind_rows(.id="cluster") %>%
  filter(grepl("Cortex L", cluster), gene %in% sig_genes) %>%
  select(cluster,gene,phi_hr) %>% 
  mutate(Age="14 months") %>% dplyr::rename(wt=phi_hr))


app=rbind(readRDS("APP_young_harmonic/par_im/coefs.rds") %>%
            bind_rows(.id="cluster") %>%
            filter(grepl("Cortex L", cluster), gene %in% sig_genes) %>%
            dplyr::select(cluster,gene,phi_hr) %>% 
            mutate(Age="7 months") %>% dplyr::rename(app=phi_hr),
          
          readRDS("APP_old_harmonic/par_im/coefs.rds") %>%
            bind_rows(.id="cluster") %>%
            filter(grepl("Cortex L", cluster), gene %in% sig_genes) %>%
            dplyr::select(cluster,gene,phi_hr) %>% 
            mutate(Age="14 months") %>%
            dplyr::rename(app=phi_hr))

coefs=left_join(wt,app, by=c("cluster", "gene", "Age"))

coefs$Age=relevel(factor(coefs$Age), "7 months")

ggplot(coefs)+geom_point(aes(x=wt, y=app, col=cluster))+
  facet_wrap(~Age)+
  geom_abline(intercept = 0,slope = 1)+
  scale_color_manual(values=cluster_color, name="")+
  xlab("NTG Phase (hrs)")+
  ylab("APP23 Phase (hrs)")+
  theme_bw(base_size = 16)+coord_fixed()+xlim(c(0,24))+ylim(c(0,24))


####################
#all
#############
setwd("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth")


wt=readRDS("WT_harmonic/par_im/coefs.rds") %>%
           bind_rows(.id="cluster") %>%
           filter(grepl("Cortex L", cluster), gene %in% sig_genes) %>%
           select(cluster,gene,phi_hr) %>% 
            dplyr::rename(wt=phi_hr)


app=readRDS("APP_harmonic/par_im/coefs.rds") %>%
            bind_rows(.id="cluster") %>%
            filter(grepl("Cortex L", cluster), gene %in% sig_genes) %>%
            dplyr::select(cluster,gene,phi_hr) %>% 
            dplyr::rename(app=phi_hr)

coefs=left_join(wt,app, by=c("cluster", "gene"))

ggplot(coefs)+geom_point(aes(x=wt, y=app, col=cluster))+
  geom_abline(intercept = 0,slope = 1)+scale_color_manual(values=cluster_color)+
  theme_bw(base_size = 16)+coord_fixed()
  
