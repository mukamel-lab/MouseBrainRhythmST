library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(glue)
library(circular)
setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(agg_exp)=names(short_names)

meta.data_org = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir))

meta.data_org$t_c=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org$t_s=as.numeric(gsub("ZT", "", meta.data_org$time))
meta.data_org= meta.data_org%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                                       t_c=round(cos(2*pi*(t_c)/24)))

rg=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/shuffle_phase/rg2.rds")

for(i in 1:50){
  
  set.seed(i)
  
  meta.data=split(meta.data_org,meta.data_org$time)
 
  meta.data=lapply(meta.data, function(df){
    set.seed(i)
    df$genotype=df$genotype[sample.int(nrow(df), replace = F)]
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
        
        maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
        gns=unique(c(counts1$gene[maxs > 3],rg[[cl]])) 
        counts2 = counts1 %>% dplyr::filter(gene %in% gns)
        
        exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
        exp_gns = exp_gns[exp_gns > 0.7]
        
        md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
        md$sex=factor(md$sex)
        md$genotype=factor(md$genotype)
        
        deseq = DESeqDataSetFromMatrix(
          countData = counts2[, c("gene", md$id)]  ,
          colData = md,
          design =  ~1,
          tidy = T
        ) %>%
          estimateSizeFactors(type = "ratio") 
        
        
        dds_7 = deseq[rg[[cl]], deseq$age=="7 months"]
        
        
        design(dds_7)= ~ sex + genotype +t_s + t_c+genotype:(t_s+t_c)
        
        dds_7 = DESeq(
          dds_7,
          test = "LRT",
          reduced = ~ sex
        )
        coefs_7 = coefficients(dds_7)
        
        coefs_7=coefs_7 %>% as.data.frame() %>% rownames_to_column("gene")%>%
          mutate( age="7 months", 
                  app_phi= atan2(t_s, t_c) %% (2 * pi),
                  wt_phi=atan2((t_s+genotypeWT.t_s), (t_c+genotypeWT.t_c)) %% (2 * pi),
                  phi_dif=atan2(sin(app_phi-wt_phi),cos(app_phi-wt_phi)),
                  phi_dif_hr=(12*phi_dif/pi)) %>% group_by(age) %>% 
          summarise(mean_phi_dif_hr=as.numeric(
            mean.circular(
              circular(phi_dif_hr, units = "hours", modulo = "2pi"),
              na.rm = TRUE
            )))
        
        
        
        dds_14 = deseq[rg[[cl]], deseq$age=="14 months"]
        
        
        design(dds_14)= ~ sex + genotype +t_s + t_c+genotype:(t_s+t_c)
        
        dds_14 = DESeq(
          dds_14,
          test = "LRT",
          reduced = ~ sex
        )
        coefs_14 = coefficients(dds_14)
        
        coefs_14=coefs_14 %>% as.data.frame() %>% rownames_to_column("gene")%>%
          mutate( age="14 months", 
                  app_phi= atan2(t_s, t_c) %% (2 * pi),
                  wt_phi=atan2((t_s+genotypeWT.t_s), (t_c+genotypeWT.t_c)) %% (2 * pi),
                  phi_dif=atan2(sin(app_phi-wt_phi),cos(app_phi-wt_phi)),
                  phi_dif_hr=(12*phi_dif/pi)) %>% group_by(age) %>% 
          summarise(mean_phi_dif_hr=as.numeric(
            mean.circular(
              circular(phi_dif_hr, units = "hours"),
              na.rm = T
            )))
        
        rbind(coefs_7, coefs_14) %>% mutate(mean_phi_dif_hr=ifelse(mean_phi_dif_hr>=12, mean_phi_dif_hr-24, mean_phi_dif_hr))
         
        
          
      })
    })
  }, mc.allow.recursive = T, ignore.interactive = T)
  
  

  res=deseqs %>% 
    bind_rows(.id="cluster") 
  
  
  
  

  saveRDS(res, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/shuffle_phase/shuffle/res{i}.rds"))
  
}
