library(tidyverse)
library(openxlsx)
library(pbmcapply)
library(DESeq2)
library(glmGamPoi)
library(customfuncs)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[, 1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

meta.data$t_c = as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s = as.numeric(gsub("ZT", "", meta.data$time))
meta.data = meta.data %>% mutate(t_s = round(sin(2 * pi * (t_s) / 24)),
                                 t_c = round(cos(2 * pi * (t_c) / 24)))
short_names2=readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(agg_exp)=plyr::mapvalues(names(agg_exp), short_names2, gsub("/","", names(short_names2)))


agg_exp = lapply(names(agg_exp), function(y) {
  x = agg_exp[[y]]
  colnames(x)[2:ncol(x)] = paste0(colnames(x)[2:ncol(x)], "__", y)
  x
})
agg_exp2 = agg_exp %>%
  purrr::reduce(function(x, y)
    left_join(x, y, by = "gene"))

md_new = data.frame(id = colnames(agg_exp2)[2:ncol(agg_exp2)]) %>%
  separate(
    id,
    into = c("sample", "region"),
    sep = "__",
    remove = F
  )

md_new = left_join(md_new, meta.data)

rhyth_g = readRDS("~/desp1/supplemental_tables/rhythmic_genes_NTG_cluster_level.rds") 
names(rhyth_g)=gsub("\\.", "", names(rhyth_g))

comparisns = combn(length(agg_exp),2)
comparisns=data.frame(Var1=names(rhyth_g)[comparisns[1,]],
                      Var2=names(rhyth_g)[comparisns[2,]])


deseqs = pbmclapply(1:nrow(comparisns), function(i) {
    try({
    r1=comparisns[i,1]
    r2=comparisns[i,2]
    
    counts1 = agg_exp2
    
    md = md_new
 
    
    md = md %>% dplyr::filter(genotype == "WT",
                     region %in% c(r1, r2))
    rg_n=c(rhyth_g[[r1]][["gene"]], rhyth_g[[r2]][["gene"]])
    counts2=counts1 %>%
      dplyr::filter(gene %in% rg_n)
    exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
    exp_gns = exp_gns[exp_gns > 0.8]
    md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
    
    deseq = DESeqDataSetFromMatrix(
      countData = counts2[, c("gene", md$id)],
      colData = md,
      design =  ~age+sex+ region + t_s + t_c + region:(t_s + t_c),
      tidy = T
    ) %>%
      estimateSizeFactors(type = "ratio") %>%
      DESeq(test = "LRT", reduced =  ~ age+sex+ region + t_s + t_c, fitType = "local")
    
    
    
    
   res=left_join(results(deseq, tidy = T) %>%
      dplyr::rename(gene=row),
  coefficients(deseq) %>%   
      as.data.frame() %>%      
      rownames_to_column("gene") ) %>% dplyr::select(-log2FoldChange,-lfcSE,-stat) %>% 
     dplyr::rename(pval_ZT_region_interaction_LRT=pvalue, FDR_BH=padj) %>% arrange(FDR_BH) %>%
     mutate(region_1=r1, region_2=r2)
  colnames(res)=gsub(r1, "R1", colnames(res))  
  colnames(res)=gsub(r2, "R2", colnames(res))  
  res[, c(1,13:14,2:12)] 
  
  })
}, ignore.interactive = T)

names(deseqs)=paste0(comparisns$Var1,"_", comparisns$Var2)


saveRDS(deseqs, "~/desp1/supplemental_tables/region_DRG_NTG_cluster_level.rds")
deseqs=deseqs[sapply(deseqs, class)=="data.frame"]
deseqs=bind_rows(deseqs)

openxlsx::write.xlsx(deseqs, "~/desp1/supplemental_tables/region_DRG_NTG_cluster_level.xlsx")
