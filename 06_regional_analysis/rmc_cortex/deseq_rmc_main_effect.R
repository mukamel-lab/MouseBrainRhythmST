library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)

setwd("~/desp1/precast/precast_final_with_ros_caud/")
data_dir = "/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("agg_c.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[, 1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

meta.data$t_c = as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s = as.numeric(gsub("ZT", "", meta.data$time))
meta.data = meta.data %>% mutate(t_s = round(sin(2 * pi * (t_s) / 24)),
                                 t_c = round(cos(2 * pi * (t_c) / 24)))
short_names =readRDS("short_names.rds")

names(agg_exp) = names(short_names)
agg_exp=agg_exp[grep("_", names(agg_exp))]
agg_exp = lapply(names(agg_exp), function(y) {
  x = agg_exp[[y]]
  colnames(x)[2:ncol(x)] = paste0(colnames(x)[2:ncol(x)], "__", y)
  x
})
agg_exp2 = agg_exp %>%
  purrr::reduce(function(x, y)
    left_join(x, y, by = "gene"))

md_new = data.frame(id = colnames(agg_exp2)[2:ncol(agg_exp2)]) %>%
  separate(
    id,
    into = c("sample", "region"),
    sep = "__",
    remove = F
  ) %>% 
  separate(region, into = c("cluster", "region"), sep = "_")

md_new = left_join(md_new, meta.data)
md_new = md_new %>% dplyr::filter(genotype=="WT")


res=pbmclapply(unique(md_new$cluster) %>% setNames(.,.), function(clst){  
  
  md = md_new %>% dplyr::filter(cluster==clst)
  counts1 = agg_exp2
  
  
  
  
  
  counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
  
  maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
  counts2 = counts1[maxs > 4, ]
  
  exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
  exp_gns = exp_gns[exp_gns > 0.8]
  
  md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
  md$region=factor(md$region)
  md$sex=factor(md$sex)
  md$age=factor(plyr::mapvalues(md$age, c("7 months", "14 months"), c("Y","O")))
  
  
  deseq = DESeqDataSetFromMatrix(
    countData = counts2[, c("gene", md$id)],
    colData = md,
    design =  ~ sex+age+region,
    tidy = T
  ) %>%
    estimateSizeFactors(type = "ratio")%>%
    DESeq(fitType = "local")
 
  
}, ignore.interactive = T)



saveRDS(res, 
        "~/desp1/precast/precast_final_with_ros_caud/analysis2/main_effect/res_region.rds",
        compress = T)
