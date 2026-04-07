library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)


setwd("~/desp1/precast/precast_final_with_ros_caud/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("agg_c.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))

agg_exp=agg_exp[grepl("^Cortex", names(agg_exp))]

deseqs=pbmclapply(agg_exp, function(cl){
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1=cl
      
      
      md = meta.data %>% dplyr::rename(id = sample) %>%
        dplyr::filter(genotype == "WT")
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      

      DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)] %>% 
                               filt(!(gene %in% c("humanAPP", "Thy1") )), 
                             colData=md, 
                             design=~sex+age+t_s+t_c,
                             tidy = T) %>%
        estimateSizeFactors(type="ratio") %>%
        estimateDispersions(fitType="local") %>%
        DESeq(test="LRT", reduced=~sex+age, fitType = "local")
      
      
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)


deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]

res=lapply(deseqs, function(x) {results(x, tidy = T)}) %>% 
  bind_rows(.id="cluster") %>%
  dplyr::rename(gene=row)

coefs=lapply(deseqs, 
             function(x) { coefficients(x) %>%
                 as.data.frame() %>% 
                 rownames_to_column("gene") 
             }) %>% 
  bind_rows(.id="cluster")
res=left_join(res,coefs)
saveRDS(res, "~/desp1/precast/precast_final_with_ros_caud/analysis2/rhythmicity/coef_wt.rds")
