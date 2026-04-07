library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)
library(Matrix.utils)

prec_run="~/desp1/precast/prec_c25q25g3000/"
setwd(prec_run)

data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

cluster_order = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")
region_order = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/region_order.rds")

all_counts= readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/fig4_app_basic/dam_scr/all_counts.rds")
seuInt=CreateSeuratObject(all_counts)
rm(all_counts)
gc()
md=  readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/fig4_app_basic/dam_scr/md_ser_all.rds") %>%
  rownames_to_column("rn")
md$anot=gsub("Olfactory Tubercle", "Amygdala", as.character(md$anot))
md = md %>% mutate(region=plyr::mapvalues(anot, cluster_order, region_order))

stopifnot(all.equal(md$rn, rownames(seuInt@meta.data)))
seuInt$region=md$region
seuInt$sample=md$sample
seuInt$region=factor(seuInt$region, levels =unique(seuInt$region)[c(1,3,2,5,4)] )
seuInt=SplitObject(seuInt, split.by = "region")
gc()
save.image("temp_reg.Rdata")
agg_c=pbmclapply(seuInt, function(smp){


  x=aggregate.Matrix(t(smp@assays$RNA@counts), groupings = smp$sample, fun="sum") %>%
    as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("gene")

  return(x)

})

saveRDS(agg_c, "objects/agg_c_regional_full.rds", compress = F)
