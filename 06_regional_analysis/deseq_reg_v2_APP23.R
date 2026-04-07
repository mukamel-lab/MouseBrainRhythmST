
#.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")

library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(glue)
setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir = "/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[, 1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

meta.data$t_c = as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s = as.numeric(gsub("ZT", "", meta.data$time))
meta.data = meta.data %>% mutate(t_s = round(sin(2 * pi * (t_s) / 24)),
                                 t_c = round(cos(2 * pi * (t_c) / 24)))
short_names=readRDS("objects/short_names2.rds")
names(short_names)[1]="L23"
agg_exp=agg_exp[short_names]
names(agg_exp)=names(short_names)

agg_exp = lapply(names(agg_exp), function(y) {
  x = agg_exp[[y]]
  colnames(x)[2:ncol(x)] = paste0(colnames(x)[2:ncol(x)], "_", y)
  x
})

agg_exp2 = agg_exp %>%
  purrr::reduce(function(x, y)
    left_join(x, y, by = "gene"))

md_new = data.frame(id = colnames(agg_exp2)[2:ncol(agg_exp2)]) %>%
  separate(
    id,
    into = c("sample", "region"),
    sep = "_",
    remove = F
  )

md_new = left_join(md_new, meta.data)

comparisns = combn(23,2)
comparisns=data.frame(Var1=names(short_names)[comparisns[1,]],
                      Var2=names(short_names)[comparisns[2,]])

regional_rhyth=pbmclapply(names(short_names) %>% setNames(.,.), function(i)
{
  counts1 = agg_exp2
  to_keep=grep(i, colnames(agg_exp2), value = T)
  counts1=counts1[,c("gene", to_keep)]
  md = md_new %>% dplyr::filter(id %in% colnames(counts1), genotype!="WT")
  counts1=counts1[,c("gene", md$id)]
  
  maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
  counts2 = counts1[maxs > 8, ]
  
  exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
  exp_gns = exp_gns[exp_gns > 0.8]
  
  md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
  
  deseq = DESeqDataSetFromMatrix(
    countData = counts2[, c("gene", md$id)]  ,
    colData = md,
    design =  ~age+sex+t_s + t_c,
    tidy = T
  ) %>%
    estimateSizeFactors(type = "ratio") %>%
    estimateDispersions(fitType = "local") %>%
    DESeq(test = "LRT",
          reduced =  ~ age+sex,
          fitType = "local")
  left_join(results(deseq, tidy = T) %>% 
    dplyr::rename(gene=row),
    coefficients(deseq) %>% as.data.frame() %>% rownames_to_column("gene") %>% 
      mutate(amp = sqrt(t_s ^ 2 + t_c ^ 2), 
             phi= atan2(t_s, t_c) %% (2 * pi),
             phi_hr=(12 / pi) * phi))

    })
  
saveRDS(regional_rhyth, "deseq_reg_int2/regional_rhythms_APP23.rds")

################################

deseqs = pbmclapply(1:nrow(comparisns) %>% setNames(., paste0(comparisns$Var1, "_", comparisns$Var2)), function(i) {
  try({
  r1=comparisns[i,1]
  r2=comparisns[i,2]
  counts1 = agg_exp2
  to_keep=grep(glue("_{r1}|_{r2}"), colnames(agg_exp2), value = T)
  counts1=counts1[,c("gene", to_keep)]
  
  md = md_new %>% dplyr::filter(id %in% colnames(counts1), genotype!="WT")
  counts1=counts1[,c("gene", md$id)]
  
  ##############################
  
  maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
  counts2 = counts1[maxs > 8, ]
  
  exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
  exp_gns = exp_gns[exp_gns > 0.8]
  
  md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
  
  deseq = DESeqDataSetFromMatrix(
    countData = counts2[, c("gene", md$id)]  ,
    colData = md,
    design =  ~age+sex+t_s + t_c,
    tidy = T
  )

  dds_r1 = deseq[, deseq$region==r1] %>%
    estimateSizeFactors(type = "ratio") %>%
    estimateDispersions(fitType = "local") %>%
    DESeq(test = "LRT",
          reduced =  ~ age+sex,
          fitType = "local")
  
  dds_r2 = deseq[, deseq$region==r2] %>%
    estimateSizeFactors(type = "ratio") %>%
    estimateDispersions(fitType = "local") %>%
    DESeq(test = "LRT",
          reduced =  ~ age+sex,
          fitType = "local")
  
  deseq = DESeqDataSetFromMatrix(
    countData = counts2[, c("gene", md$id)]  ,
    colData = md %>% mutate(region=factor(region, levels=c(r1,r2))),
    design =  ~ age+sex+region  + t_s + t_c+region:(t_s+t_c),
    tidy = T
  ) %>%
    estimateSizeFactors(type = "ratio") %>%
    estimateDispersions(fitType = "local") %>%
    DESeq(test = "LRT",
          reduced =  ~ age+ sex + region  + t_s + t_c,
          fitType = "local")
  r1g=results(dds_r1, tidy = T) %>% dplyr::filter(padj<0.1)%>% .[["row"]]
  r2g=results(dds_r2, tidy = T) %>% dplyr::filter(padj<0.1)%>% .[["row"]]
  res=results(deseq, tidy=T) %>%
    dplyr::filter(row %in% c(r1g, r2g))%>%
    dplyr::rename(gene=row) %>%
    mutate(rhythmic_in_r1=(gene %in% r1g),
           rhythmic_in_r2=(gene %in% r2g))
  
  res$fdr=p.adjust(res$pvalue, method = "BH")
  
  coefs = coefficients(deseq) %>% as.data.frame() %>% rownames_to_column("gene")
  
  colnames(coefs)=gsub(r1, "r1", colnames(coefs))
  colnames(coefs)=gsub(r2, "r2", colnames(coefs))
  
  coefs = coefs %>% mutate(r1_amp = sqrt(t_s ^ 2 + t_c ^ 2), 
                           r2_amp=sqrt((t_s+regionr2.t_s) ^ 2 + (t_c+regionr2.t_c) ^ 2),
                           r1_phi= atan2(t_s, t_c) %% (2 * pi),
                           r2_phi=atan2(t_s+ regionr2.t_s, t_c+ regionr2.t_c) %% (2 * pi),
                           r1_phi_hr=(12 / pi) * r1_phi,
                           r2_phi_hr=(12 / pi) * r2_phi)
                        

  coefs$region_1=r1
  coefs$region_2=r2
  
  
  res = left_join(res, coefs) 
  r1g=results(dds_r1, tidy = T) %>% 
    dplyr::select(row,padj) %>%
    dplyr::rename(gene=row, rhythmicity_padj_r1=padj)
  r2g=results(dds_r2, tidy = T) %>% 
    dplyr::select(row,padj) %>% 
    dplyr::rename(gene=row, rhythmicity_padj_r2=padj)
  
  rg12=left_join(r1g,r2g)
  res=left_join(res, rg12)%>%
    dplyr::relocate(
      gene,
      region_1,
      region_2,
      fdr,
      rhythmic_in_r1,
      rhythmic_in_r2,
      rhythmicity_padj_r1,
      rhythmicity_padj_r2,
      r1_amp,
      r2_amp,
      r1_phi,
      r2_phi,
      r1_phi_hr,
      r2_phi_hr
    )
  
  
  return(res)
  })
}, ignore.interactive = T)

saveRDS(deseqs,"deseq_reg_int2/regional_DRG_APP23_v2.rds")
