library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(openxlsx)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir = "/cndd2/agelber/hal/qc_aligned"


agg_exp = readRDS("objects/agg_c.rds")
names(agg_exp) = gsub(' ', "_", names(agg_exp))
names(agg_exp) = make.names(names(agg_exp))
nm_bu=names(agg_exp)

agg_exp = lapply(names(agg_exp), function(nm) {
  df = agg_exp[[nm]]
  
  colnames(df)[2:ncol(df)] = paste0(colnames(df)[2:ncol(df)], "_", nm)
  
  df
})

names(agg_exp)=nm_bu

agg_exp_full = agg_exp[[1]]

for (i in 2:length(agg_exp)) {
  agg_exp_full = left_join(agg_exp_full, agg_exp[[i]])
}

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[, 1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

meta.data$t_c = as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s = as.numeric(gsub("ZT", "", meta.data$time))
meta.data = meta.data %>% mutate(t_s = round(sin(2 * pi * (t_s) / 24)),
                                 t_c = round(cos(2 * pi * (t_c) / 24)))


meta.data2 = data.frame(sample = colnames(agg_exp_full)[2:ncol(agg_exp_full)])
meta.data2 = meta.data2 %>% separate(sample,
                                     c("sample", "region"),
                                     sep = "_",
                                     extra = "merge")
meta.data = left_join(meta.data2, meta.data)


counts1 = agg_exp_full

md = meta.data %>% dplyr::mutate(id = paste0(sample, "_", region)) %>%
  relocate(id, .before = everything())


maxs = apply(counts1 %>% dplyr::select(-gene) , 1, max)
counts2 = counts1[maxs > 10,]
cs = colSums(counts2 %>% select(-gene))
cs = names(cs)[cs > 2e5]


ngens = apply(counts2 %>% dplyr::select(-gene), 2, function(x) {
  sum(x > 0)
})
ngens = ngens / nrow(counts2)
ngens = names(ngens)[ngens > 0.9]

md$sex = relevel(factor(md$sex), "M")
md$genotype = relevel(factor(md$genotype), "WT")
md$age = relevel(factor(md$age), "7 months")

md = md %>% filt(id %in% intersect(cs, ngens), genotype=="WT",
                 region %in% names(agg_exp)[c(1:5,11)])

r_genes = left_join(
  readRDS(
    "~/desp1/precast/prec_c25q25g3000/deseq_rhyth2/WT_harmonic/par_im/res.rds"
  ),
  readRDS(
    "~/desp1/precast/prec_c25q25g3000/deseq_rhyth2/WT_harmonic/par_im/coefs.rds"
  ) %>%
    bind_rows(.id = "cluster")
) %>% mutate(rel_amp = amp / log2(baseMean + 1)) %>%
  dplyr::filter(cluster %in% c("Cortex Layer 2/3", "Dentate Gyrus-sg"),
                padj < 0.05,
                rel_amp > 0.03)
md$region[grepl("Cortex", md$region)]="Cortex"
deseq = DESeqDataSetFromMatrix(
  countData = counts2[, c("gene", md$id)] %>%
    filt(!(gene %in% c("humanAPP", "Thy1"))) %>%
    dplyr::filter(gene %in% r_genes$gene),
  colData = md,
  design =  ~age+sex+t_s+t_c+region:(t_s+t_c),
  tidy = T
) %>%
  estimateSizeFactors(type = "poscounts") %>%
  estimateDispersions(fitType = "glmGamPoi") %>%
  DESeq(test="LRT", reduced=~age+sex+t_s+t_c)




res = results(deseq,
          tidy = T) %>%
  dplyr::rename(gene = row)

res_sig = res %>% dplyr::filter(padj < 0.1)
