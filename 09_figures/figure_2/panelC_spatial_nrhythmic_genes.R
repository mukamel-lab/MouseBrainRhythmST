library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)
library(scattermore)

prec_run="~/desp1/precast/prec_c17q17g3000/"
setwd(prec_run)

data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

seuInt=readRDS("seuInt.rds")
seuInt$pos1=seuInt@reductions$position@cell.embeddings[,1]
seuInt$pos2=seuInt@reductions$position@cell.embeddings[,2]

md=seuInt@meta.data

md=md %>% dplyr::filter(sample=="107-B")
res_sig_all=readRDS("deseq/res_sig_app.rds")
res_sig_all = res_sig_all %>% group_by(cluster) %>% tally
res_sig_all=res_sig_all %>%dplyr::rename(anot=cluster)
md=left_join(md, res_sig_all, by="anot")
md=md %>% na.replace(0)
p=ggplot(md)+ geom_point(aes(x=pos1,y=pos2, col=n))+theme_void()+
  ggtitle(glue("Number of DEG Per Cluster"))+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_fixed() + scale_x_reverse()+scale_y_reverse()
p

cl_size=seuInt@meta.data %>% group_by(anot) %>% tally %>% dplyr::rename(cl_size=n)
md=left_join(md, cl_size, by="anot")

p1=ggplot(md)+ geom_point(aes(x=pos1,y=pos2, col=n/cl_size))+theme_void()+
  ggtitle(glue("Number of DEG Per Cluster Scaled by Cluster Size"))+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_fixed() + scale_x_reverse()+scale_y_reverse()
p1
p|p1
