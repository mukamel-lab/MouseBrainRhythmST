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
  dplyr::filter(sample %in% list.files(data_dir)) %>% dplyr::filter(age=="7 months")

meta.data_org$t_c=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org$t_s=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org= meta.data_org%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                                       t_c=round(cos(2*pi*(t_c)/24)))


  meta.data=meta.data_org
  
  deseqs=pbmclapply(agg_exp, function(cl){
    # full model specified in design
    # reduced model specified in reduced
    try({
      suppressWarnings({
        counts1=cl
        
        md = meta.data %>% dplyr::rename(id = sample)
        
        counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
        
        maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
        counts2 = counts1[maxs > 3, ]
        
        exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
        exp_gns = exp_gns[exp_gns > 0.7]
        
        md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
        md$sex=factor(md$sex)
        md$genotype=factor(md$genotype)
        
        deseq = DESeqDataSetFromMatrix(
          countData = counts2[, c("gene", md$id)]  ,
          colData = md,
          design =  ~ genotype+sex + t_s + t_c+genotype:(t_s+t_c),
          tidy = T
        ) %>%
          estimateSizeFactors(type = "ratio") 
        rhyth_g=results(deseq %>% 
                          DESeq(
                            test = "LRT",
                            reduced = ~ sex+genotype
                          ), tidy=T) %>% 
          dplyr::filter(padj<0.05)%>%
          pull(row)
        #######
        dds_w = deseq[, deseq$genotype == "WT"]
        dds_w$genotype=droplevels(dds_w$genotype)
        design(dds_w)= ~ sex + t_s + t_c
        
        dds_w = DESeq(
          dds_w,
          test = "LRT",
          reduced = ~ sex
        )
        coefs_w = coefficients(dds_w)
        
        coefs_w=coefs_w[rhyth_g, ] %>% as.data.frame() %>% rownames_to_column("gene")%>%
          mutate(amp = sqrt(t_s ^ 2 + t_c ^ 2), genotype="WT")
        ###########################################
        dds_a = deseq[, deseq$genotype == "APP23"]
        dds_a$genotype=droplevels(dds_a$genotype)
        design(dds_a)= ~ sex + t_s + t_c
        
        dds_a = DESeq(
          dds_a,
          test = "LRT",
          reduced = ~ sex
        )
        coefs_a = coefficients(dds_a) 
        coefs_a=coefs_a[rhyth_g, ] %>% as.data.frame()  %>% rownames_to_column("gene")%>%
          mutate(amp = sqrt(t_s ^ 2 + t_c ^ 2),genotype="APP23") 
        
        
        list(res=
               left_join(coefs_a, coefs_w, by="gene", suffix = c(".a",".w")) %>% mutate(dif=amp.a-amp.w) %>%
               summarise(mean_dif=mean(dif)),
             rg=rhyth_g)
       
        
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  
  res=lapply(deseqs, function(y) {y$res}) %>% 
    bind_rows(.id="cluster") 
  
saveRDS(res, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/shuffle_amp/true_young.rds"))
  

rg=lapply(deseqs, function(y) {y$rg}) 

saveRDS(rg, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/shuffle_amp/rg_young.rds"))
