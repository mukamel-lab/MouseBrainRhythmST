library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))

res_f_m = pbmclapply(names(agg_exp), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample) %>%
        dplyr::filter(genotype == "WT")
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "poscounts")
      
      dds_f = deseq[, deseq$sex == "F"]
      
      dds_f = DESeq(
        dds_f,
        test = "LRT",
        reduced = ~ 1
      )
     
      coefs_f=left_join(coefficients(dds_f) %>% 
                          as.data.frame() %>% 
                          rownames_to_column("gene"),
                        results(dds_f, tidy = T) %>% dplyr::rename(gene = row) %>%
                          mutate(sex ="F"))
      counts_f=
      dds_m = deseq[, deseq$sex == "M"]
      
      dds_m = DESeq(
        dds_m,
        test = "LRT",
        reduced = ~ 1
      )
      
      
      coefs_m=left_join(coefficients(dds_m) %>% 
                          as.data.frame() %>% 
                          rownames_to_column("gene"),
                        results(dds_m, tidy = T) %>% dplyr::rename(gene = row) %>%
                          mutate(sex ="M"))
      
      ncounts_m = counts(dds_m, normalized = T) %>% as.data.frame() %>%
        rownames_to_column("gene")  %>%
        pivot_longer(cols = c(-gene),
                     names_to = "id",
                     values_to = "l2expr") %>%
        mutate(l2expr = log2(l2expr + 1)) %>% left_join(md)
      
      ncounts_f = counts(dds_f, normalized = T) %>% as.data.frame() %>%
        rownames_to_column("gene")  %>%
        pivot_longer(cols = c(-gene),
                     names_to = "id",
                     values_to = "l2expr") %>%
        mutate(l2expr = log2(l2expr + 1)) %>% left_join(md)
      
      list(coefs = list("F"=coefs_f, M=coefs_m), counts = list("F"=ncounts_f, M=ncounts_m))
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)

names(res_f_m) = names(agg_exp)

saveRDS(res_f_m,"~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/for_plotting/counts_coefs_plot.rds")

