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



deseqs = pbmclapply(names(agg_exp), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample) %>%
        dplyr::filter(age == "7 months")
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ genotype + sex + t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio") %>%
        estimateDispersions(fitType = "parametric") %>%
        DESeq(test = "LRT",
              reduced =  ~ genotype+sex,
              fitType = "parametric")
      
      res = results(deseq, tidy = T) %>%
        dplyr::rename(gene = row)
      
      
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ genotype + sex + t_s + t_c + genotype:(t_s + t_c),
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio") %>%
        estimateDispersions(fitType = "local") %>%
        DESeq(test = "LRT",
              reduced =  ~ genotype + sex + t_s + t_c,
              fitType = "local")
      
 
      
      coefs = coefficients(deseq) %>% as.data.frame() %>% rownames_to_column("gene")
      coefs = coefs %>% mutate(app_amp = sqrt(t_s ^ 2 + t_c ^ 2), 
                               app_phi= atan2(t_s, t_c) %% (2 * pi),          
                               app_phi_hr=(12 / pi) * app_phi,
                               wt_amp = sqrt((t_s + genotypeWT.t_s) ^ 2 + (t_c + genotypeWT.t_c) ^ 2),
                               
                               wt_phi=atan2(t_s+ genotypeWT.t_s, t_c+ genotypeWT.t_c) %% (2 * pi),
                               wt_phi_hr=(12 / pi) * wt_phi
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
res_a_w = pbmclapply(names(agg_exp), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample) %>%
        dplyr::filter(age == "7 months")
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
      deseq = DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ sex + t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio")
      
      dds_a = deseq[, deseq$genotype == "APP23"]
      
      dds_a = DESeq(
        dds_a,
        test = "LRT",
        reduced = ~ sex,
        fitType = "local"
      )
      res_a = results(dds_a, tidy = T) %>% dplyr::rename(gene = row) %>%
        mutate(genotype ="APP23")
      
      dds_w = deseq[, deseq$genotype == "WT"]
      
      dds_w= DESeq(
        dds_w,
        test = "LRT",
        reduced = ~ sex,
        fitType = "local"
      )
      res_w = results(dds_w, tidy = T) %>% dplyr::rename(gene = row) %>% 
        mutate(genotype= "WT")
      res_aw = rbind(res_w, res_a)
      res_aw
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)

names(res_a_w) = names(agg_exp)

saveRDS(coefs, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/coefs_y.rds")

saveRDS(ncounts, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/counts_y.rds")

saveRDS(res_a_w %>% bind_rows(.id="cluster"),"~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/young/res_y_aw.rds")

