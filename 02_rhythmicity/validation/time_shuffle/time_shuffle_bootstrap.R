library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)


setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c_ds_25spots.rds")

meta_orig = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 




bstrp=1:30
names(bstrp)=paste0("s", bstrp)
# bst_res=lapply(bstrp, function (i) {
#   
#   meta.data=meta_orig
#   
#   meta.data$time=meta.data$time[sample.int(nrow(meta.data), nrow(meta.data), replace = F)]
#   
#   
#   meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
#   meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
#   meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
#                                  t_c=round(cos(2*pi*(t_c)/24)))
# deseqs=pbmclapply(agg_exp, function(cl){
#   # full model specified in design
#   # reduced model specified in reduced
#   try({
#     suppressWarnings({
#       counts1=cl
#       
#       md=meta.data %>% dplyr::rename(id=sample)
#       maxs=apply(counts1 %>% dplyr::select(-gene) , 1, max)
#       counts2=counts1[maxs>10, ]
#       cs=colSums(counts2 %>% dplyr::select(-gene))
#       cs=names(cs)[cs>2e5]
#       ngens=apply(counts2 %>% dplyr::select(-gene), 2, function(x){sum(x>0)})
#       ngens=ngens/nrow(counts2)
#       ngens=names(ngens)[ngens>0.9]
#       
#       md=md %>% filt(id %in% intersect(cs, ngens), genotype!="WT")
#       
#       DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)] %>% 
#                                filt(!(gene %in% c("humanAPP", "Thy1") )), 
#                              colData=md, 
#                              design=~age+sex+t_s+t_c,
#                              tidy = T) %>%
#         estimateSizeFactors(type="poscounts") %>%
#         estimateDispersions(fitType="glmGamPoi") %>%
#         DESeq(test="LRT", reduced=~age+sex)
#       
#       
#     })
#   })
# }, mc.allow.recursive = T, ignore.interactive = T)
# 
# 
# deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]
# 
# res=lapply(deseqs, function(x) {results(x, tidy = T)}) %>% 
#   bind_rows(.id="cluster") %>%
#   dplyr::rename(gene=row)
# res_sig= res %>% dplyr::filter(padj<0.05)
# res_sig %>% group_by(cluster) %>% tally
# })
# 
# app_res=bst_res %>% bind_rows(.id="run")

###########################################
bst_res=lapply(bstrp, function (i) {
  
  meta.data=meta_orig
  
  meta.data$time=meta.data$time[sample.int(nrow(meta.data), nrow(meta.data), replace = F)]
  
  
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
        cs=colSums(counts2 %>% dplyr::select(-gene))
        cs=names(cs)[cs>2e5]
        ngens=apply(counts2 %>% dplyr::select(-gene), 2, function(x){sum(x>0)})
        ngens=ngens/nrow(counts2)
        ngens=names(ngens)[ngens>0.9]
        
        md=md %>% filt(id %in% intersect(cs, ngens), genotype=="WT")
        
        DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)] %>% 
                                 filt(!(gene %in% c("humanAPP", "Thy1") )), 
                               colData=md, 
                               design=~age+sex+t_s+t_c,
                               tidy = T) %>%
          estimateSizeFactors(type="poscounts") %>%
          estimateDispersions(fitType="glmGamPoi") %>%
          DESeq(test="LRT", reduced=~age+sex)
        
        
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  
  deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]
  
  res=lapply(deseqs, function(x) {results(x, tidy = T)}) %>% 
    bind_rows(.id="cluster") %>%
    dplyr::rename(gene=row)
  res_sig= res %>% dplyr::filter(padj<0.05)
  res_sig %>% group_by(cluster) %>% tally
})

wt_res=bst_res %>% bind_rows(.id="run")


setwd("deseq_rhth/time_shuffle/")
#dir.create("bootstrap")
setwd("bootstrap")

saveRDS(wt_res, "wt_res.rds")
#saveRDS(app_res, "app_res.rds")