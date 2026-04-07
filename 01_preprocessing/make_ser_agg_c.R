library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)
library(Matrix.utils)
library(Matrix)

prec_run="~/desp1/precast/prec_c25q25g3000/"
setwd(prec_run)

data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

seuInt=readRDS("seuInt_qc.rds")

seuInt=SplitObject(seuInt, "sample")

ser_objs=pbmclapply(names(seuInt), function(s) {
  x=Load10X_Spatial(glue("{data_dir}/{s}/outs")) %>%
    subset(cells=seuInt[[s]]$cell)
  x$sample=s
  md=x@meta.data
  md$cell=rownames(md)
  md=left_join(md, seuInt[[s]]@meta.data)
  md=left_join(md, metadata)
  md=as.data.frame(md)
  rownames(md)=md$cell
  
  
  x=CreateSeuratObject(x@assays$Spatial@counts)
  x@meta.data=md
  x=SCTransform(x)
  return(x)
  
})
names(ser_objs) = metadata$sample

var_feat=unique(unlist(lapply(ser_objs, VariableFeatures)))

ser_objs=merge(ser_objs[[1]], ser_objs[2:65])



VariableFeatures(ser_objs)=var_feat


ser_objs_split=SplitObject(ser_objs, "anot")

agg_c=pbmclapply(ser_objs_split, function(smp){


  x=aggregate.Matrix(t(smp@assays$RNA@counts), groupings = smp$sample, fun="sum") %>%
    as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("gene")

  return(x)

})
saveRDS(agg_c, "objects/agg_c.rds", compress = F)

ser_objs=RunPCA(ser_objs)

ser_objs=RunUMAP(ser_objs, dims=1:30)
ser_objs=FindNeighbors(ser_objs, dims=1:30, reduction="pca")
ser_objs=FindClusters(ser_objs)

saveRDS(ser_objs, "ser_all.rds")
saveRDS(ser_objs@meta.data, "md_ser_all.rds")

DefaultAssay(ser_objs)="RNA"
ser_objs=NormalizeData(ser_objs) %>% ScaleData() %>% FindVariableFeatures()

Idents="anot"

marks=FindAllMarkers(ser_objs, logfc.threshold = .5)

saveRDS(marks, "marks.rds")

