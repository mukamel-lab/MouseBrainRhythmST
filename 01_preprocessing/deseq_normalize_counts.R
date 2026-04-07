library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(openxlsx)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

  agg_exp = readRDS("objects/agg_c.rds")
  
meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
    mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
    dplyr::filter(sample %in% list.files(data_dir)) 
  
meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


deseqs=pbmclapply(agg_exp, function(cl){
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1=cl
    
      md=meta.data %>% dplyr::rename(id=sample)
      maxs=apply(counts1 %>% dplyr::select(-gene) , 1, max)
      counts2=counts1[maxs>10, ]
      cs=colSums(counts2 %>% select(-gene))
      cs=names(cs)[cs>2e5]
      ngens=apply(counts2 %>% dplyr::select(-gene), 2, function(x){sum(x>0)})
      ngens=ngens/nrow(counts2)
      ngens=names(ngens)[ngens>0.9]
      
      md$sex=relevel(factor(md$sex), "M")
      md$genotype=relevel(factor(md$genotype), "WT")
      md$age=relevel(factor(md$age), "7 months")
      
      md=md %>% filt(id %in% intersect(cs, ngens))
      
      DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)], 
                             colData=md, 
                             design=~sex+age,
                             tidy = T) %>%
        estimateSizeFactors(type="poscounts")
      
      
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)


deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]

ncounts=lapply(deseqs, 
               function(x) {counts(x, normalized=T) %>%
                   as.data.frame()})

saveRDS(ncounts, "objects/agg_counts_deseq_norm.rds")

nlcpm=lapply(ncounts, 
               function(x) {x1=x+1; log2(x1)})

saveRDS(nlcpm, "objects/deseq_norm_l2cpm.rds")

