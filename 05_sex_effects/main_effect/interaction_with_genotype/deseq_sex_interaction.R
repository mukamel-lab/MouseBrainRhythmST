library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")
# agg_exp = readRDS("objects/agg_c_ds_25spots.rds")
cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


rgs=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/main_effect/interaction_with_genotype/res_main_by_sex_young.rds")

rgs=rgs %>% dplyr::filter(padj<0.1)
rgs=split(rgs, rgs$cluster)

rgs=lapply(rgs, function(y) {unique(y$gene)})
deseqs = pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample) %>% dplyr::filter(age=="7 months")
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ sex +genotype+sex:genotype,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio") 
      
      deseq[rgs[[cl_n]], ]
      
     
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)





saveRDS(deseqs, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/main_effect/interaction_with_genotype/deseq_interaction_young.rds")
