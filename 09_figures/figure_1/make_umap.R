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

setwd("~/desp1/precast/prec_c25q25g3000/figures/figure_1_plots/")

data_dir="/cndd2/agelber/hal/qc_aligned"
metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

seuInt=readRDS("~/desp1/precast/prec_c25q25g3000/seuInt_qc.rds")
seuInt$anot=factor(as.character(seuInt$anot),
                   levels = readRDS("../cluster_order.rds"))
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


################################
# umaps
################################
dir.create("umaps")
setwd("umaps")

cols_cluster = c(ggsci::pal_material("green")(7)[2:7],
                 ggsci::pal_material("blue")(7)[2:7],
                 ggsci::pal_material("orange")(4)[2:4],
                 
                 ggsci::pal_material("red")(6)[2:6],
                 ggsci::pal_material("deep-purple")(4)[2:4]
                 )
names(cols_cluster)=levels(seuInt$anot)

p1=ggplot(md)+ 
  geom_scattermore(interpolate = T, 
    aes(x=umap1,y=umap2, col=anot))+
  theme_bw()+
  guides(col=guide_legend(
    title = "Annotated Cluster", ncol=2,
    override.aes = list(size=4)))+
  ggtitle("UMAP on PRECAST Embedding")+
  scale_color_manual(values=cols_cluster)+
  geom_text(data=md %>% group_by(anot) %>% summarise(umap1=mean(umap1), umap2=mean(umap2)),
            aes(x=umap1, y=umap2,label=anot))+
  xlab("UMAP1")+ylab("UMAP2")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        panel.grid = element_blank(),
        axis.text = element_blank()
          )


spatplts=lapply(names(md_spl), function(s){
  df=md_spl[[s]]
  ggplot(df)+ geom_point(aes(y=-1*imagerow,x=imagecol, col=anot))+
    theme_void()+
    guides(col="none")+
    ggtitle(glue("Slice {s}"))+
    scale_color_manual(values=cols_cluster)+
    theme(axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.title = element_blank(),
          panel.grid = element_blank(),
          title = element_text(size=18))+
    coord_fixed()
  
  
})

names(spatplts)=names(md_spl)
p=spatplts$`030-B`
pan_a=(p1/p)+plot_layout(guides="collect")

pdf("umap_slice.pdf", width = 20, height=12)
pan_a
dev.off()


sex=ggplot(md %>% mutate(sex=plyr::mapvalues(sex, c("F","M"), c("Female", "Male"))))+ 
  geom_scattermore(interpolate = T, 
                   aes(x=umap1,y=umap2, col=sex))+
  theme_bw()+
  guides(col=guide_legend(
    title = "", ncol=1,
    override.aes = list(size=4)))+
  ggtitle("Sex")+
  ggsci::scale_color_npg()+
  xlab("")+ylab("")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+coord_fixed()


age=ggplot(md)+ 
  geom_scattermore(interpolate = T, 
                   aes(x=umap1,y=umap2, col=age))+
  theme_bw()+
  guides(col=guide_legend(
    title = "", ncol=1,
    override.aes = list(size=4)))+
  ggtitle("Age")+
  scale_color_manual(values=ggsci::pal_npg()(5)[c(3,5)])+
  xlab("")+ylab("")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+coord_fixed()

geno=ggplot(md %>% mutate(genotype=plyr::mapvalues(genotype, c("APP23", "WT"), c("APP23-TG", "NTG"))))+ 
  geom_scattermore(interpolate = T, 
                   aes(x=umap1,y=umap2, col=genotype))+
  theme_bw()+
  guides(col=guide_legend(
    title = "", ncol=1,
    override.aes = list(size=4)))+
  scale_color_manual(values=ggsci::pal_jco()(7)[c(6,7)])+
  ggtitle("Genotype")+
  xlab("")+ylab("")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+coord_fixed()


time=ggplot(md)+ 
  geom_scattermore(interpolate = T, 
                   aes(x=umap1,y=umap2, col=time))+
  theme_bw()+
  guides(col=guide_legend(
    title = "", ncol=1,
    override.aes = list(size=4)))+
  ggtitle("Zeitgeber Time")+
  ggsci::scale_color_igv()+
  xlab("")+ylab("")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16))+coord_fixed()

pan_cov=wrap_elements(((age|sex)/(time | geno)))+
  plot_annotation(title = "UMAPs Colored by Covariates")& 
  theme(plot.title = element_text(hjust=.5 ,size=20, face="bold"))
  
pdf("umap_covariates.pdf", width = 10, height = 8)
pan_cov
dev.off()

