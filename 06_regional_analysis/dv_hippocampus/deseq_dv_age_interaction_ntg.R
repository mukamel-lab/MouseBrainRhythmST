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
clusts=str_split_fixed(names(agg_exp), "_",2)[,1] %>% unique()
names(clusts)=clusts

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

deg = readRDS("analysis/deseq_rhyth/res_NTG.rds") %>% 
  dplyr::filter(padj<0.1)
deg$cluster=plyr::mapvalues(deg$cluster, short_names, names(short_names))
deg=split(deg,deg$age)

deg = lapply(deg, function(df) {split(df, df$cluster)})
md_new=md_new %>% dplyr::filter(genotype=="WT")
md_new=split(md_new, md_new$age)

res=lapply(md_new, function(mdta){
  
 deseqs=pbmclapply(clusts, function(i) {
  try({
    
    smps=intersect(grep(i, colnames(agg_exp2), value=T), mdta$id)
    
    counts1 = agg_exp2[,c('gene', smps)]
    
    md = mdta %>% dplyr::filter(id %in% colnames(counts1))
    
    degs=deg[[md$age[1]]][grep(i, names(deg[[md$age[1]]]), value=T)] %>% bind_rows()
    
    counts2=counts1 %>% dplyr::filter(gene %in% degs$gene)
    
    md=md %>% mutate(region=str_split_fixed(region, "_",2)[,2]) %>%
      mutate(region=factor(region, levels=c("V","D")))
    
    exp_gns=colSums(counts2 %>% dplyr::select(-gene)>0)/nrow(counts2)
    
    exp_gns=exp_gns[exp_gns>=0.8]
    
    md=md %>% dplyr::filter(id %in% intersect(colnames(counts2), names(exp_gns)))
    
    counts2=counts2[,c("gene",intersect(md$id, colnames(counts2)))]
    
     DESeqDataSetFromMatrix(
      countData = as.data.frame(counts2) ,
      colData = as.data.frame(md),
      design =  ~ sex+region+t_s+t_c+region:(t_s+t_c),
      tidy = T
    ) %>%
      estimateSizeFactors(type = "poscounts") %>%
      estimateDispersions(fitType = "glmGamPoi") %>%
      DESeq(tes="LRT", reduced = ~ sex+region+t_s+t_c)
    
  })
}, ignore.interactive = T)
 
 deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]
 
 
 lapply(deseqs, function(x) {results(x, tidy = T) }) %>% 
   bind_rows(.id="cluster") %>%
   dplyr::rename(gene=row)
}) %>% bind_rows(.id="age")

saveRDS(res, "analysis/deseq_rhyth_dv_interact/hip_NTG_reg_int.rds")


