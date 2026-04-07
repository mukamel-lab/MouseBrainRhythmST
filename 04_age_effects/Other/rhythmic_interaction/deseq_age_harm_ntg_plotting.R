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

res_o_y = pbmclapply(names(agg_exp), function(cl_n) {
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
      
      dds_o = deseq[, deseq$age == "14 months"]
      
      dds_o = DESeq(
        dds_o,
        test = "LRT",
        reduced = ~ 1
      )
     
      coefs_o=left_join(coefficients(dds_o) %>% 
                          as.data.frame() %>% 
                          rownames_to_column("gene"),
                        results(dds_o, tidy = T) %>% dplyr::rename(gene = row) %>%
                          mutate(age ="14 months"))

      dds_y = deseq[, deseq$age == "7 months"]
      
      dds_y = DESeq(
        dds_y,
        test = "LRT",
        reduced = ~ 1
      )
      
      
      coefs_y=left_join(coefficients(dds_y) %>% 
                          as.data.frame() %>% 
                          rownames_to_column("gene"),
                        results(dds_y, tidy = T) %>% dplyr::rename(gene = row) %>%
                          mutate(age ="7 months"))
      
      ncounts_o = counts(dds_o, normalized = T) %>% as.data.frame() %>%
        rownames_to_column("gene")  %>%
        pivot_longer(cols = c(-gene),
                     names_to = "id",
                     values_to = "l2expr") %>%
        mutate(l2expr = log2(l2expr + 1)) %>% left_join(md)
      
      ncounts_y = counts(dds_y, normalized = T) %>% as.data.frame() %>%
        rownames_to_column("gene")  %>%
        pivot_longer(cols = c(-gene),
                     names_to = "id",
                     values_to = "l2expr") %>%
        mutate(l2expr = log2(l2expr + 1)) %>% left_join(md)
      
      list(coefs = list(O=coefs_o, Y=coefs_y), counts = list(O=ncounts_o, Y=ncounts_y))
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)

names(res_o_y) = names(agg_exp)

saveRDS(res_o_y,"~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/age_int_new/for_plotting/counts_coefs_plot.rds")

