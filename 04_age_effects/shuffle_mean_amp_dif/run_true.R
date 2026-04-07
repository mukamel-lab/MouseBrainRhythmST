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
# meta.data_org = meta.data_org %>%
#   group_by(sex, time) %>%
#   group_modify(~ {
#     m <- min(table(.x$age))
#     .x %>%
#       group_by(age, .add = TRUE) %>%
#       slice_sample(n = m)
#   }) %>%
#   ungroup()

  deseqs=pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl){
    # full model specified in design
    # reduced model specified in reduced
    try({
      suppressWarnings({
        counts1=agg_exp[[cl]]
        
        md = meta.data_org %>% dplyr::rename(id = sample) 
        
        counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
        
        maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
        counts2 = counts1[maxs >4, ]
        
        exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
        exp_gns = exp_gns[exp_gns > 0.8]
        
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

        rgs=deseq %>% 
          DESeq(
                test = "LRT",
                reduced = ~ sex,
                fitType = "parametric"
          ) %>% results(tidy=T) %>% 
             slice_min(order_by = padj, n = 200) %>% pull(row)
        
        dds_y=deseq[rgs,deseq$age=="7 months"] 
        
        dds_y= DESeq(dds_y,
            test = "LRT",
            reduced = ~ sex,
            fitType = "parametric"
          )
        

        
        dds_o=deseq[rgs,deseq$age=="14 months"]
        dds_o=dds_o%>% 
          DESeq(
            test = "LRT",
            reduced = ~ sex,
            fitType = "parametric"
          )
        
        
        coefs_y=coefficients(dds_y) %>% as.data.frame() %>% rownames_to_column("gene")  %>%
          mutate(amp=sqrt(t_s ^ 2 + t_c ^ 2), age="Y") %>%
          group_by(age) %>% summarise(amp=mean(amp))
        
        coefs_o=coefficients(dds_o) %>% as.data.frame() %>% rownames_to_column("gene") %>%
           mutate(amp=sqrt(t_s ^ 2 + t_c ^ 2), age="O")%>%
          group_by(age) %>% summarise(amp=mean(amp))
        
        coefs=rbind(coefs_o,coefs_y)
        list(coefs=coefs, rgs=rgs)
        
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  
  
  res=lapply(deseqs, function(x) {x$coefs}) %>% 
    bind_rows(.id="cluster")

  saveRDS(deseqs, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
                    "age_int_new/shuffle_mean_amp_dif/rg_OY.rds"))
  
  saveRDS(res, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
                    "age_int_new/shuffle_mean_amp_dif/true_OY.rds"))

