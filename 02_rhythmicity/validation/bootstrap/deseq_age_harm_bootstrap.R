#!/usr/bin/env Rscript
chnks = commandArgs(trailingOnly=TRUE)

library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(glue)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_app/ds_bootstrap/")

data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp_list = readRDS("agg_c_ds_10spots_20min_100samples.rds")
setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth/ds_bootstrap/")
dir.create("WT_age")
setwd("WT_age")

dir.create("14_months_objects")
setwd("14_months_objects")

names(agg_exp_list)=paste0("smp_", 1:length(agg_exp_list))

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


for (smp in names(agg_exp_list)) {
  
  try({
    deseqs=pbmclapply(agg_exp_list[[smp]], function(cl){
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
          
          md=md %>% filt(id %in% intersect(cs, ngens), age=="14 months", genotype=="WT")
          
          DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)] %>% 
                                   filt(!(gene %in% c("humanAPP", "Thy1") )), 
                                 colData=md, 
                                 design=~sex+t_c+t_s,
                                 tidy = T) %>%
            estimateSizeFactors(type="poscounts") %>%
            estimateDispersions(fitType="glmGamPoi") %>%
            DESeq(test="LRT", reduced=~sex)
          
          
        })
      })
    }, mc.allow.recursive = T, ignore.interactive = T)
    
    
    deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]
    
    res=lapply(deseqs, function(x) {results(x, tidy = T)}) %>% 
      bind_rows(.id="cluster") %>%
      dplyr::rename(gene=row)
    res_sig= res %>% subset(padj<0.05) 
    saveRDS(res_sig, paste0(smp, ".rds"))
  })
}

setwd("..")

dir.create("7_months_objects")
setwd("7_months_objects")



for (smp in names(agg_exp_list)) {
  
  try({
    deseqs=pbmclapply(agg_exp_list[[smp]], function(cl){
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
          
          md=md %>% filt(id %in% intersect(cs, ngens), age=="7 months", genotype=="WT")
          
          DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)] %>% 
                                   filt(!(gene %in% c("humanAPP", "Thy1") )), 
                                 colData=md, 
                                 design=~sex+t_s+t_c,
                                 tidy = T) %>%
            estimateSizeFactors(type="poscounts") %>%
            estimateDispersions(fitType="glmGamPoi") %>%
            DESeq(test="LRT", reduced=~sex)
          
          
        })
      })
    }, mc.allow.recursive = T, ignore.interactive = T)
    
    
    deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]
    
    res=lapply(deseqs, function(x) {results(x, tidy = T)}) %>% 
      bind_rows(.id="cluster") %>%
      dplyr::rename(gene=row)
    res_sig= res %>% subset(padj<0.05) 
    saveRDS(res_sig, paste0(smp, ".rds"))
  })
}
