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
meta.data_org <- meta.data_org %>% 
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
        
        md = meta.data_org %>% dplyr::rename(id = sample) 
        
        counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
        
        maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
        counts2 = counts1[maxs >4, ]
        
        exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
        exp_gns = exp_gns[exp_gns > 0.8]
        
        md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
       
        
        deseq = DESeqDataSetFromMatrix(
          countData = counts2[, c("gene", md$id)]  ,
          colData = md,
          design =  ~age+t_s+t_c,
          tidy = T
        ) %>%
          estimateSizeFactors(type = "ratio") 
        
      
        
        dds_f = deseq[, deseq$sex == "F"]
        
        dds_f = DESeq(
          dds_f,
          test = "LRT",
          reduced = ~ age,
          fitType = "parametric")
        
        
        # res_f=results(dds_f, tidy=T)  %>%   slice_min(order_by = padj, n = 100)
        # 
        # coefs_f=coefficients(dds_f) %>% as.data.frame() %>% rownames_to_column("gene") %>%
        #     dplyr::filter(gene %in% res_f$row) %>% mutate(amp=sqrt(t_s ^ 2 + t_c ^ 2), sex="F") %>% 
        #   group_by(sex) %>% summarise(amp=mean(amp))
        res_f=results(dds_f, tidy=T)  %>%
          summarise(n1=sum(padj<0.1, na.rm = T),
                    n05=sum(padj<0.05, na.rm = T),
                    n01=sum(padj<0.01, na.rm = T))
        dds_m=deseq[,deseq$sex=="M"] %>% 
          DESeq(
            test = "LRT",
            reduced = ~ age,
            fitType = "parametric"
          )
        
        # res_m=results(dds_m, tidy=T) %>%   slice_min(order_by = padj, n = 100)
        # coefs_m=coefficients(dds_m) %>% as.data.frame() %>% rownames_to_column("gene") %>%
        #   dplyr::filter(gene %in% res_m$row) %>% mutate(amp=sqrt(t_s ^ 2 + t_c ^ 2), sex="M") %>% 
        #   group_by(sex) %>% summarise(amp=mean(amp))
        # coefs=rbind(coefs_f, coefs_m)  
        res_m=results(dds_m, tidy=T) %>%
        summarise(n1=sum(padj<0.1, na.rm = T),
                    n05=sum(padj<0.05, na.rm = T),
                    n01=sum(padj<0.01, na.rm = T))
        
        res=list(Male=res_m, Female=res_f) %>% bind_rows(.id="sex")
        
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  
  
  res=deseqs %>% 
    bind_rows(.id="cluster") 
  
  
  saveRDS(res, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
                    "sex_int_geno_spec/shuffle_num_rg/true_MF_balanced.rds"))

