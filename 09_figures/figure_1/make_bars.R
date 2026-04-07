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

bar_by_samp=ggplot(md 
                   , aes(x=anot, fill=sample, col=genotype))+
  geom_bar(position = "fill")+ 
  scale_y_continuous(labels = scales::percent)+
  coord_flip()+
  scale_fill_manual(name="Sample", values = c(ggsci::pal_igv()(51), 
                                              ggsci::pal_d3()(10),
                                              ggsci::pal_material()(4)))+
  scale_color_manual(name="Genotype", values = ggsci::pal_d3()(2))+
  theme(text = element_text(size=23))+theme_bw()+
  xlab("Cluster")+ylab("Percentage")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text= element_text(size=15))+
  ggtitle("Percentage of Sample Per Cluster")



# pdf("bar_by_samp.pdf", width = 15, height = 9)
# bar_by_samp
# dev.off()

#############################################
# bar by clust
#############################################


bar_by_clust=ggplot(md, aes(x=sample, fill=anot))+
  geom_bar(position = "fill")+ 
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(name="Cluster", values = cols_cluster)+
  scale_color_manual(name="Genotype", 
                     values = ggsci::pal_d3()(2))+
  theme(text = element_text(size=23))+theme_bw()+
  xlab("Cluster")+
  facet_grid(~ genotype+age+sex,
             scales = "free_x",
             space = "free_x",
             switch = "x") +
  theme(title = element_text(size=18), 
        axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text.y= element_text(size=15),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, units = "cm"))+
  ggtitle("Percent Cluster by Sample")+
  ylab("Percentage")
  
# pdf("bar_by_clust.pdf", width = 15, height = 9)
# bar_by_clust
# dev.off()

#############################################
# spots per clust
#############################################


spots_per_clust=md %>% group_by(genotype, anot) %>% tally()
spc_ord=spots_per_clust %>% dplyr::filter(genotype=="NTG") %>% arrange(n)
spots_per_clust$anot=factor(spots_per_clust$anot, levels = spc_ord$anot)


spc_plt=ggplot(spots_per_clust  %>% 
                 mutate(anot=factor(anot, levels = levels(seuInt$anot))))+geom_col(aes(x=anot, fill=genotype, y=n), 
                                                                                   position = position_dodge())+
  ggsci::scale_fill_d3(name="Genotype")+theme_bw()+
  xlab("Cluster")+ylab("Number of Spots")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text.y= element_text(size=15),
        axis.text.x= element_text(size=15,angle=70,hjust=1 ))+
  ggtitle("Spots per Cluster")

# pdf("spots_per_clust.pdf", width = 15, height = 9)
# spc_plt
# dev.off()

#############################################
# combine
#############################################


bar_by_samp1=bar_by_samp+
  theme(plot.margin =
          margin(t = 2, r = 0,
                 b = 2,
                 l = 2,
                 unit = "pt"))+
  guides(fill="none",
         color=guide_legend(
           override.aes=list(lwd=7,fill=NA)))

spc_plt1=spc_plt+coord_flip()+
  theme(axis.text.y = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 0, unit = "pt"))+
  xlab("")+guides(fill="none")

pct_spts_per=ggpubr::ggarrange(bar_by_samp1, 
                       spc_plt1, 
                       ncol = 2, nrow = 1,
                       align = "h", 
                       common.legend = T,
                       legend = "right")

pct_spts_per
