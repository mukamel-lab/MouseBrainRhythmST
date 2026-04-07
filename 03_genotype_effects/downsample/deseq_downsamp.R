library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)
library(glue)
setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


for(i in 1:30){
  agg_exp=readRDS(glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
                       "geno_int_new/downsample/spot10/agg_c_ds{i}.rds"))

deseqs = pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample) 
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 5)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.7]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      md$sex=factor(md$sex)
      md$age=factor(plyr::mapvalues(md$age, c("7 months", "14 months"), c("Y","O")))
      md$genotype=factor(md$genotype)

      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ sex + t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio")
    
      dds_ao = deseq[, deseq$genotype == "APP23" & deseq$age=="O"]
      
      res_ao = results(
        dds_ao %>% estimateDispersions(fitType = "local") %>%
          DESeq(
            test = "LRT",
            reduced = ~ sex,
            fitType = "local"
          ), tidy=T
      ) %>% dplyr::rename(gene = row) %>%
        mutate(genotype ="APP23", age="O")
     
      dds_wy = deseq[, deseq$genotype == "WT"  & deseq$age == "Y"]
      res_wy = results(
        dds_wy %>% estimateDispersions(fitType = "local") %>%
          DESeq(
            test    = "LRT",
            reduced = ~ sex,
            fitType = "local"
          ),
        tidy = TRUE
      ) %>%
        dplyr::rename(gene = row) %>%
        mutate(genotype = "WT", age = "Y")
      
      
      # Genotype = "WT", Age = "O"
      dds_wo = deseq[, deseq$genotype == "WT"  & deseq$age == "O"]
      res_wo = results(
        dds_wo %>% estimateDispersions(fitType = "local") %>%
          DESeq(
            test    = "LRT",
            reduced = ~ sex,
            fitType = "local"
          ),
        tidy = TRUE
      ) %>%
        dplyr::rename(gene = row) %>%
        mutate(genotype = "WT", age = "O")
      
      
      # Genotype = "APP23", Age = "Y"
      dds_ay = deseq[, deseq$genotype == "APP23" & deseq$age == "Y"]
      res_ay = results(
        dds_ay %>% estimateDispersions(fitType = "local") %>%
          DESeq(
            test    = "LRT",
            reduced = ~ sex,
            fitType = "local"
          ),
        tidy = TRUE
      ) %>%
        dplyr::rename(gene = row) %>%
        mutate(genotype = "APP23", age = "Y")
      
    bind_rows(list(res_ao, res_ay, res_wo, res_wy))  %>% group_by(age, genotype) %>%
      summarise(n0.1=sum(padj<0.1,na.rm=T),n0.01=sum(padj<0.01,na.rm=T), n0.05=sum(padj<0.05,na.rm=T))
    
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T) %>% bind_rows(.id="cluster")


saveRDS(deseqs, glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
                     "geno_int_new/downsample/spot10_res/res{i}.rds"))


}