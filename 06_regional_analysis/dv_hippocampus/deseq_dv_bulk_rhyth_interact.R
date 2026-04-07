library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)

setwd("~/desp1/precast/precast_final_with_RMC_DV/")
data_dir = "/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("agg_c.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[, 1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

meta.data$t_c = as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s = as.numeric(gsub("ZT", "", meta.data$time))
meta.data = meta.data %>% mutate(t_s = round(sin(2 * pi * (t_s) / 24)),
                                 t_c = round(cos(2 * pi * (t_c) / 24)))
short_names =readRDS("short_names.rds")

agg_exp=agg_exp[grep("Dorsal|Ventral", names(agg_exp), value = T)]
sn=names(short_names)
names(sn)=short_names
names(agg_exp) =as.character(sn[names(agg_exp)])

agg_exp = lapply(names(agg_exp), function(y) {
  x = agg_exp[[y]]
  reg=str_split_fixed(y, "_",2)[,2]
  clust=str_split_fixed(y, "_",2)[,1]
  colnames(x)[2:ncol(x)] = paste0(reg,"__",colnames(x)[2:ncol(x)], "___", clust)
  x
})
agg_exp2 = agg_exp %>%
  purrr::reduce(function(x, y)
    left_join(x, y, by = "gene"))


agg_exp3 = agg_exp2 %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>% 
  mutate(prefix = gsub("___.*","", sample)) %>%
  group_by(gene, prefix) %>%
  dplyr::summarise(expr = sum(expr), .groups = "drop")  %>%
  mutate(prefix = gsub("__","_", prefix)) %>%
  pivot_wider(id_cols = gene, names_from = prefix, values_from = expr)


md_new = data.frame(id = colnames(agg_exp3)[2:ncol(agg_exp3)]) %>%
  separate(
    id,
    into = c("region", "sample"),
    sep = "_",
    remove = F
  )

md_new = left_join(md_new, meta.data)
md_new=md_new %>% dplyr::filter(genotype=="WT")


rhyth_g = readRDS("analysis/deseq_rhyth/hip_regional/res_NTG.rds") %>% 
  dplyr::filter(padj<0.1) 

rhyth_g = split(rhyth_g, rhyth_g$age)
md_new=split(md_new, md_new$age)


res=lapply(md_new, function(mdta) {
  

      md=mdta
      counts1 = agg_exp3[,c("gene",intersect(md$id, colnames(agg_exp3)))]
      
      rt_gns=rhyth_g[[mdta$age[1]]] 
      rt_gns=unique(rt_gns$gene)
      
      counts2=counts1 %>% dplyr::filter(gene %in% rt_gns)
      
      exp_gns=colSums(counts2 %>% dplyr::select(-gene)>0)/nrow(counts2)
      exp_gns=exp_gns[exp_gns>=0.8]
      md=md %>% dplyr::filter(id %in% intersect(colnames(counts2), names(exp_gns)))
      
      
      counts2=counts2[,c("gene",intersect(md$id, colnames(counts2)))]
      
      deseq = DESeqDataSetFromMatrix(
        countData = as.data.frame(counts2), 
        colData = as.data.frame(md),
        design =  ~ region + t_s + t_c + region:(t_s + t_c),
        tidy = T
      ) %>%
        estimateSizeFactors(type = "poscounts") %>%
        estimateDispersions(fitType = "glmGamPoi") %>%
        DESeq(test = "LRT", reduced =  ~ region + t_s + t_c)
      
      
      
      
      left_join(results(deseq, tidy = T) %>%
                      dplyr::rename(gene=row),
                    coefficients(deseq) %>%   
                      as.data.frame() %>%      
                      rownames_to_column("gene"))
      
     
})

res=res %>% bind_rows(.id="age")

saveRDS(res, "analysis/deseq_rhyth_dv_interact/hip_regional_NTG_reg_int.rds")
