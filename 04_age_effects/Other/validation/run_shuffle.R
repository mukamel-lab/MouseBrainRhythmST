.libPaths("~/R/x86_64-pc-linux-gnu-library/4.3/")
library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(glue)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")


meta.data_org = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) %>% dplyr::filter(genotype=="WT")

meta.data_org$t_c=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org$t_s=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org= meta.data_org%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                                       t_c=round(cos(2*pi*(t_c)/24)))


run_one_shuffle= function(i) {
  
  
  set.seed(i)
  
  balanced <- meta.data_org %>% 
    group_by(sex, time) %>%              # work one block at a time
    group_modify(~ {                     # .x = data, .y = keys for this block
      m <- min(table(.x$age))            # smallest of the 2 sex-specific counts
      .x %>% 
        group_by(age, .add = TRUE) %>%   # now sample *within* each sex
        slice_sample(n = m)              # keeps m rows at random
    }) %>% 
    ungroup()
  

  
  meta.data=split(balanced,paste0(balanced$sex,balanced$time))
 
  meta.data=lapply(meta.data, function(df){
    set.seed(i)
    df$age=df$age[sample.int(nrow(df), replace = F, )]
    df
  })
  meta.data=bind_rows(meta.data)
  
  deseqs=pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl){
    # full model specified in design
    # reduced model specified in reduced
    try({
      suppressWarnings({
        counts1=agg_exp[[cl]]
        
        md = meta.data %>% dplyr::rename(id = sample) 
        
        counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
        
        maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 1)
        counts2 = counts1[maxs >=4, ]
        
        exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
        
        exp_gns_cut=sample(c(0.7,0.75, 0.8), 1)
        exp_gns = exp_gns[exp_gns > exp_gns_cut]
        
        md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
 
        md$sex=factor(md$sex)
        md$genotype=factor(md$genotype)
        
        deseq = DESeqDataSetFromMatrix(
          countData = counts2[, c("gene", md$id)]  ,
          colData = md,
          design =  ~sex+t_c+t_s,
          tidy = T
        ) %>%
          estimateSizeFactors(type = "ratio") 
        dds_y=deseq[,deseq$age=="7 months"] %>% 
          DESeq(
            test = "LRT",
            reduced = ~ sex,
            fitType = "parametric"
          )
        
        res_y=results(dds_y, tidy=T)%>%  
          summarise(n1=sum(padj<0.1, na.rm = T),
                    n05=sum(padj<0.05, na.rm = T),
                    n01=sum(padj<0.01, na.rm = T))
        
        dds_o=deseq[,deseq$age=="14 months"] %>% 
          DESeq(
            test = "LRT",
            reduced = ~ sex,
            fitType = "parametric"
          )
        
        res_o=results(dds_o, tidy=T) %>%  
          summarise(n1=sum(padj<0.1, na.rm = T),
                    n05=sum(padj<0.05, na.rm = T),
                    n01=sum(padj<0.01, na.rm = T))
        
        res=list(O=res_o, Y=res_y) %>% bind_rows(.id="age")
     
     
     })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  

  res=deseqs %>% 
    bind_rows(.id="cluster")  %>%
    pivot_wider(id_cols = c( "cluster"), names_from = "age", 
                                           values_from = starts_with("n")) %>%
    mutate(n1_rat=log2(n1_O/n1_Y),
           n05_rat=log2(n05_O/n05_Y), 
           n01_rat=log2(n01_O/n01_Y))
  

  

  saveRDS(res, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
  "age_int_new/shuffle_num_rg/shuffles/res{i}.rds"))
  invisible(NULL)
}
pbmclapply(1:100, run_one_shuffle, mc.cores = 6, mc.allow.recursive = T)