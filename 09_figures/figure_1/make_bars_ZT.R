library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)
library(pbmcapply)
library(scattermore)
library(egg)

setwd("~/desp1/precast/prec_c25q25g3000/figures/figure_1_plots/")

data_dir="/cndd2/agelber/hal/qc_aligned"
metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 


seuInt=readRDS("~/desp1/precast/prec_c25q25g3000/seuInt_qc.rds")


seuInt$anot=factor(as.character(seuInt$anot),
                   levels = rev(readRDS("../cluster_order.rds")))
md=seuInt@meta.data
md$umap1=seuInt@reductions$umap@cell.embeddings[,1]
md$umap2=seuInt@reductions$umap@cell.embeddings[,2]

md_spl=split(md, md$sample)

md=pbmclapply(names(md_spl), function(s) {
  x=Load10X_Spatial(glue("{data_dir}/{s}/outs")) %>%
    subset(cells=md_spl[[s]]$cell)
  y=x@meta.data
  y$cell=rownames(y)
  y=cbind(y, x@images$slice1@coordinates)
  
  
  
  y=left_join(md_spl[[s]], y)
  return(y)
  
})

names(md)=names(md_spl)
md=bind_rows(md, .id="sample")

md=left_join(md, metadata)
md_spl=split(md, md$sample)


####################
#bar plts
#######################


dir.create("bar_plt")
setwd("bar_plt/")

md_ord=md %>% select(sample, genotype) %>%
  distinct() %>%arrange(genotype)
md$sample=factor(md$sample, levels = md_ord$sample)
md$genotype=plyr::mapvalues(md$genotype,
                            c("APP23", "WT"), 
                            c("APP23-TG", "NTG"))
md$age=factor(md$age, levels = rev(unique(md$age)))
md$sex=factor(md$sex)

cols_cluster = c(ggsci::pal_material("green")(7)[2:7],
                 ggsci::pal_material("blue")(7)[2:7],
                 ggsci::pal_material("orange")(4)[2:4],
                 
                 ggsci::pal_material("red")(6)[2:6],
                 ggsci::pal_material("deep-purple")(4)[2:4]
)

names(cols_cluster)=levels(seuInt$anot)

#############################################
# bar by samp
#############################################

bar_by_samp=ggplot(md %>% 
                     mutate(Time=factor(time, levels=paste0("ZT", 6*0:3)),
                            genotype=factor(genotype, levels=rev(unique(genotype)))) 
                   , aes(x=anot, fill=Time))+
  geom_bar(position = "fill")+ 
  scale_y_continuous(labels = scales::percent)+
   coord_flip()+
  theme(text = element_text(size=23))+theme_bw()+
  xlab("Cluster")+ylab("Percentage")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text= element_text(size=15))+
  ggtitle("ZT Distribution by Cluster")+facet_wrap(~genotype)+
  ggsci::scale_fill_npg()

bar_by_samp

# pdf("bar_by_samp.pdf", width = 15, height = 9)
# bar_by_samp
# dev.off()

