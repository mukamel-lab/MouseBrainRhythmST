library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c_ds_25spots.rds")
cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))



deseqs = pbmclapply(names(agg_exp), function(cl_n) {
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
        design =  ~ age + sex + t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio") %>%
        estimateDispersions(fitType = "local") %>%
        DESeq(test = "LRT",
              reduced =  ~ age+sex,
              fitType = "local")
      
      res = results(deseq, tidy = T) %>%
        dplyr::rename(gene = row)
      
      
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ age + sex + t_s + t_c + sex:(t_s + t_c),
        tidy = T
      ) %>%
        estimateSizeFactors(type = "poscounts") %>%
        estimateDispersions(fitType = "local") %>%
        
        DESeq(test = "LRT",
              reduced =  ~ age + sex + t_s + t_c,
              fitType = "local")
      
      
      
      coefs = coefficients(deseq) %>% as.data.frame() %>% rownames_to_column("gene")
      coefs = coefs %>% mutate(f_amp = sqrt(t_s ^ 2 + t_c ^ 2), 
                               m_amp = sqrt((t_s + sexM.t_s) ^ 2 + (t_c + sexM.t_c) ^ 2),
                               f_phi= atan2(t_s, t_c) %% (2 * pi),
                               m_phi=atan2(t_s+ sexM.t_s, t_c+ sexM.t_c) %% (2 * pi),
                               f_phi_hr=(12 / pi) * f_phi,
                               m_phi_hr=(12 / pi) * m_phi
      )
      coefs = left_join(coefs, res)
      res = results(deseq, tidy = T) %>%
        dplyr::rename(gene = row)
      coefs = left_join(
        coefs,
        res %>% dplyr::select(gene, pvalue, padj) %>% dplyr::rename(pvalue_int =
                                                                      pvalue, padj_int = padj)
      )
      ncounts = counts(deseq, normalized = T) %>% as.data.frame() %>%
        rownames_to_column("gene")  %>%
        pivot_longer(cols = c(-gene),
                     names_to = "id",
                     values_to = "l2expr") %>%
        mutate(l2expr = log2(l2expr + 1)) %>% left_join(md)
      list(coefs = coefs, counts = ncounts)
      
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)

names(deseqs) = names(agg_exp)
coefs = lapply(deseqs, function(x) {x$coefs}) %>% bind_rows(.id = "cluster")
ncounts=lapply(deseqs, function(x) {x$counts}) %>% bind_rows(.id = "cluster")
res_f_m = pbmclapply(names(agg_exp), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample) %>%
        dplyr::filter(genotype == "WT")
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 5)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ age + t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio")
      
      dds_f = deseq[, deseq$sex == "F"]
      
      dds_f = DESeq(
        dds_f,
        test = "LRT",
        reduced = ~ age,
        fitType = "local",
        useT = TRUE
      )
      res_f = results(dds_f, tidy = T) %>% dplyr::rename(gene = row) %>%
        mutate(sex ="F")
      
      dds_m = deseq[, deseq$sex == "M"]
      
      dds_m = DESeq(
        dds_m,
        test = "LRT",
        reduced = ~ age,
        fitType = "local",
        useT = TRUE
      )
      res_m = results(dds_m, tidy = T) %>% dplyr::rename(gene = row) %>% mutate(sex ="M")
      res_mf = rbind(res_m, res_f)
      res_mf
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)

names(res_f_m) = names(agg_exp)

saveRDS(coefs, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/downsample/coefs_wt_ds.rds")

saveRDS(ncounts, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/downsample/counts_wt_ds.rds")

saveRDS(res_f_m %>% bind_rows(.id="cluster"),"~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/downsample/res_wt_mf_ds.rds")

