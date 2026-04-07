
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")

library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)

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
cluster_ord = readRDS("objects/cluster_order.rds")
short_names = names(agg_exp)
names(short_names) = c(
  "CTX23",
  "CTX4",
  "CTX5a",
  "CTX5b",
  "CTX6a",
  "CTXsp",
  "CA1",
  "CA3so",
  "CA3sp",
  "DGmo",
  "DGsg",
  "RHP",
  "COAa",
  "COAp",
  "Piriform",
  "Caudoputamen",
  "STRv",
  "Amygdala",
  "GP",
  "RN",
  "LV",
  "FP",
  "Meninges"
)
names(agg_exp) = names(short_names)
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
md_new = md_new %>% mutate(reg_c = ifelse(
  region %in% c("CTX23", "CTX4", "CTX5a", "CTX5b", "CTX6a",
                "CTXsp"),
  "CTX",
  region
))
rhyth_g = readRDS("deseq_rhyth2/WT_harmonic/par_im/res.rds") %>% 
  dplyr::filter(padj<0.1)

rhyth_g = split(rhyth_g, rhyth_g$cluster)

names(rhyth_g) = plyr::mapvalues(names(rhyth_g), short_names, names(short_names))
comparisns = combn(23,2)
comparisns=data.frame(Var1=names(rhyth_g)[comparisns[1,]],
                      Var2=names(rhyth_g)[comparisns[2,]])


deseqs = pbmclapply(1:nrow(comparisns), function(i) {
  try({
  r1=comparisns[i,1]
  r2=comparisns[i,2]
  counts1 = agg_exp2
  
  md = md_new
  maxs = apply(counts1 %>% dplyr::select(-gene) , 1, max)
  counts2 = counts1[maxs > 10, ]
  cs = colSums(counts2 %>% dplyr::select(-gene))
  cs = names(cs)[cs > 2e5]
  ngens = apply(counts2 %>% dplyr::select(-gene), 2, function(x) {
    sum(x > 0)
  })
  ngens = ngens / nrow(counts2)
  ngens = names(ngens)[ngens > 0.9]
  
  md = md %>% filt(id %in% intersect(cs, ngens),
                   genotype != "WT",
                   region %in% c(r1, r2))
  rg_n=c(rhyth_g[[r1]][["gene"]], rhyth_g[[r2]][["gene"]])
  
  deseq = DESeqDataSetFromMatrix(
    countData = counts2[, c("gene", md$id)] %>%
      filt(!(gene %in% c(
        "humanAPP", "Thy1"
      ))) %>%
      dplyr::filter(gene %in% rg_n),
    colData = md,
    design =  ~ region + t_s + t_c + region:(t_s + t_c),
    tidy = T
  ) %>%
    estimateSizeFactors(type = "poscounts") %>%
    estimateDispersions(fitType = "glmGamPoi") %>%
    DESeq(test = "LRT", reduced =  ~ region + t_s + t_c)
  
  
  

  res=results(deseq, tidy = T) %>%
    dplyr::rename(gene=row)
  coefs=coefficients(deseq) %>%   
    as.data.frame() %>%      
    rownames_to_column("gene") 
  return(list(coefs=coefs, res=res))
  
  })
}, ignore.interactive = T)

names(deseqs)=paste0(comparisns$Var1,"_", comparisns$Var2)

saveRDS(deseqs, 
        "~/desp1/precast/prec_c25q25g3000/deseq_wt_reg_int/all_deseqs_app.rds",
        compress = T)
