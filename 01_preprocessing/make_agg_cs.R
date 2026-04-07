library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)
library(Matrix.utils)
library(DESeq2)
library(glmGamPoi)
library(customfuncs)

prec_run="~/desp1/precast/prec_c25q25g3000/"
setwd(prec_run)

data_dir="/cndd2/agelber/hal/qc_aligned"


seuInt=readRDS("objects/ser_all.rds")

setwd("rhth_ds_bootstrap/filt50_samp02_resamp_100bs/")


data_dir = "/cndd2/agelber/hal/qc_aligned"


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[, 1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

meta.data$t_c = as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s = as.numeric(gsub("ZT", "", meta.data$time))
meta.data = meta.data %>% mutate(t_s = round(sin(2 * pi * (t_s) / 24)),
                                 t_c = round(cos(2 * pi * (t_c) / 24)))


seuInt=SplitObject(seuInt, split.by = "anot")

for(i in 1:100){
  set.seed(i)
agg_c=pbmclapply(seuInt, function(smp){
  set.seed(i)
  md=smp@meta.data %>% rownames_to_column("rn")
 
  md_min=md %>% group_by(sample) %>%  dplyr::filter(n()>=50) %>% tally
  md_min=round(min(md_min$n))
  md=md %>% group_by(sample) %>% dplyr::filter(n()>=50) %>%
    sample_n(md_min, replace = T)
  

  x=aggregate.Matrix(t(smp@assays$RNA@counts[,md$rn]), groupings = md$sample, fun="sum") %>%
    as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("gene")
  
  return(x)
  
})
saveRDS(agg_c, paste0("agg_c_", i,".rds"))
}
