library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)
library(openxlsx)

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
      
      coefs = coefficients(deseq) %>% as.data.frame() %>% rownames_to_column("gene") %>%
        mutate(amp = sqrt(t_s ^ 2 + t_c ^ 2), 
                               phi= atan2(t_s, t_c) %% (2 * pi),
                               phi_hr=(12 / pi) * phi  )
      
      left_join(res, coefs, by="gene") %>% dplyr::filter(padj<0.1)
     
      
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)

names(deseqs) = names(agg_exp)
short_names2 <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")

names(deseqs)=plyr::mapvalues(names(deseqs), short_names2, names(short_names2))

deseqs=lapply(deseqs, function(x) { x %>% dplyr::select(-log2FoldChange,-lfcSE,-stat, -phi) %>% 
    dplyr::rename(pval_rhythmicity_LRT=pvalue, FDR_BH=padj, log2_amplitude=amp, acrophase=phi_hr) %>% arrange(FDR_BH)
  
  })
names(deseqs)=make.names(names(deseqs))
saveRDS(deseqs, "~/desp1/supplemental_tables/rhythmic_genes_NTG_cluster_level.rds")
openxlsx::write.xlsx(deseqs, "~/desp1/supplemental_tables/rhythmic_genes_NTG_cluster_level.xlsx")
