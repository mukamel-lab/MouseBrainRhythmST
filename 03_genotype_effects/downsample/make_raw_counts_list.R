library(Seurat)
library(tidyverse)
library(customfuncs)
library(PRECAST)
library(glue)
library(pbmcapply)
library(Matrix.utils)
library(customfuncs)



seuInt=readRDS("~/desp1/precast/prec_c25q25g3000/objects/ser_all.rds")

precast_metadata = readRDS("~/desp1/precast/prec_c25q25g3000/objects/md_ser_all_new.rds")

md=seuInt@meta.data
md$anot=NULL
md=left_join(md %>% rownames_to_column("rn"), 
             precast_metadata %>% dplyr::select(sample,cell,anot))%>%
  column_to_rownames("rn")
md=md[rownames(seuInt@meta.data),]
seuInt$anot=md$anot



seuInt=SplitObject(seuInt, split.by = "anot")

raw_counts_list=lapply(seuInt, function(seu) {
  
  list(md=seu@meta.data, counts=seu[["RNA"]]$counts)
})

saveRDS(raw_counts_list, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/downsample/raw_counts_list.rds")
