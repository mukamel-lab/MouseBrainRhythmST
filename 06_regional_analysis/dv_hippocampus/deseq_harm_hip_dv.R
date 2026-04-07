library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)


setwd("~/desp1/precast/precast_final_with_RMC_DV/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("agg_c.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


short_names =readRDS("short_names.rds")

agg_exp=agg_exp[grep("Ventral|Dorsal", names(agg_exp), value = T)]
sn=names(short_names)
names(sn)=short_names
names(agg_exp) =as.character(sn[names(agg_exp)])

meta.data=meta.data %>% dplyr::filter(genotype=="WT")

meta.data=split(meta.data, paste0(meta.data$age))




res=lapply(meta.data, function(mdta){
  
  deseqs=pbmclapply(agg_exp, function(cl){
    # full model specified in design
    # reduced model specified in reduced
    try({
      suppressWarnings({
        
        md=mdta %>% dplyr::rename(id=sample)
        
        counts1=cl[, c("gene",intersect(md$id, colnames(cl)))]
        
        maxs=rowSums((counts1 %>% dplyr::select(-gene))>10)
        
        counts2=counts1[maxs>2, ]
        
        exp_gns=colSums(counts2 %>% dplyr::select(-gene)>0)/nrow(counts2)
        
        exp_gns=exp_gns[exp_gns>=0.8]
        
        md=md %>% dplyr::filter(id %in% intersect(colnames(counts2), names(exp_gns)))
        
        counts2=counts2[,c("gene",intersect(md$id, colnames(counts2)))]
        
        DESeqDataSetFromMatrix(countData=as.data.frame(counts2) , 
                               colData=as.data.frame(md), 
                               design=~sex+t_c+t_s,
                               tidy = T) %>%
          estimateSizeFactors(type="poscounts") %>%
          estimateDispersions(fitType="glmGamPoi") %>%
          DESeq(test="LRT", reduced=~sex)
        
        
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  
  deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]
  
  res=lapply(deseqs, function(x) {results(x, tidy = T) }) %>% 
    bind_rows(.id="cluster") %>%
    dplyr::rename(gene=row)
  
  
  
  
  coefs=lapply(deseqs, 
               function(x) { coefficients(x) %>%
                   as.data.frame() %>% 
                   rownames_to_column("gene") 
               })%>% 
    bind_rows(.id="cluster")
  
  
  
  
  res=left_join(res, coefs)
  
  
  ncounts = lapply(deseqs,
                   function(x) {
                     counts(x, normalized = T) %>%
                       as.data.frame() %>%
                       rownames_to_column("gene") %>%
                       pivot_longer(-gene, names_to = "sample", values_to = "norm_counts") %>%
                       left_join(mdta %>%
                                   dplyr::select(-t_c,-t_s)) %>%
                       mutate(l2expr = log2(norm_counts + 1))
                   }) %>%  bind_rows(.id = "cluster")
  list(res=res, counts=ncounts)
  
  
})

coefs=list("7 months" =res$`7 months`$res,
           "14 months" =res$`14 months`$res) %>% bind_rows(.id="age")

saveRDS(coefs,"analysis/deseq_rhyth/res_NTG.rds")


ncounts=list("7 months" =res$`7 months`$counts,
           "14 months" =res$`14 months`$counts) %>% bind_rows(.id="age")

saveRDS(ncounts,"analysis/deseq_rhyth/ncounts_NTG.rds")
