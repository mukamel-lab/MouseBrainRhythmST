library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(glue)
library(qvalue)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")


meta.data_org = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) %>% dplyr::filter(genotype=="WT", age=="7 months")

meta.data_org$t_c=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org$t_s=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org= meta.data_org%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                                       t_c=round(cos(2*pi*(t_c)/24)))

run_one_shuffle= function(i) {
  set.seed(i)
  
  
  balanced <- meta.data_org %>% 
    group_by(age, time) %>%              # work one block at a time
    group_modify(~ {                     # .x = data, .y = keys for this block
      m <- min(table(.x$sex))            # smallest of the 2 sex-specific counts
      .x %>% 
        group_by(sex, .add = TRUE) %>%   # now sample *within* each sex
        slice_sample(n = m)              # keeps m rows at random
    }) %>% 
    ungroup()
  
  
  
  

  
  deseqs=pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl){
    # full model specified in design
    # reduced model specified in reduced
    try({
      suppressWarnings({
        counts1=agg_exp[[cl]]
        
        md = balanced %>% dplyr::rename(id = sample) 
        
        counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
        
        maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
        counts2 = counts1[maxs >=3, ]
        
        exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
        
        exp_gns_cut=sample(c(0.7,0.75, 0.8), 1)
        exp_gns = exp_gns[exp_gns > exp_gns_cut]
        
        md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
        
        md$sex=factor(md$sex)
        md$genotype=factor(md$genotype)
        
        deseq = DESeqDataSetFromMatrix(
          countData = counts2[, c("gene", md$id)]  ,
          colData = md,
          design =  ~1+t_c+t_s,
          tidy = T
        ) %>%
          estimateSizeFactors(type = "ratio") 
        
        dds_f=deseq[,deseq$sex=="F"] 
        dds_f=dds_f%>% 
          DESeq(
            test = "LRT",
            reduced = ~ 1,
            fitType = "glmGamPoi"
          )
        
        res_f=results(dds_f, tidy=T)
        q_f=qvalue(res_f$pvalue,  pi0.method="bootstrap", lambda = seq(0,0.95,0.05))
        q_f=q_f$pi0
        q_fs=qvalue(res_f$pvalue, pi0.method="smoother", smooth.log.pi0=T, lambda = seq(0,0.95,0.05))
        q_fs=q_fs$pi0
        
        dds_m=deseq[,deseq$sex=="M"] 
        dds_m=dds_m%>% 
          DESeq(
            test = "LRT",
            reduced = ~ 1,
            fitType = "glmGamPoi"
          )
        
        res_m=results(dds_m, tidy=T) 
        
        q_m=qvalue(res_m$pvalue, pi0.method="bootstrap", lambda = seq(0,0.95,0.05))
        q_m=q_m$pi0
        q_ms=qvalue(res_m$pvalue, pi0.method="smoother", smooth.log.pi0=T, lambda = seq(0,0.95,0.05))
        q_ms=q_ms$pi0
        data.frame(seed=i, pi0_m=q_m, pi0_f=q_f, pi0_m_smooth=q_ms, pi0_f_smooth=q_fs)  
        
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  
  
  res=deseqs %>% 
    bind_rows(.id="cluster") 
  
  
  
  
  saveRDS(res, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
                    "sex_int_geno_spec/shuffle_q_value/true_young/res{i}.rds"))
  invisible(NULL)
}
pbmclapply(1:15, run_one_shuffle, mc.cores = 4, mc.allow.recursive = T)