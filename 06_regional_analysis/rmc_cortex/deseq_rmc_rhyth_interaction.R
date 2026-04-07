library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)

setwd("~/desp1/precast/precast_final_with_ros_caud/")
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

names(agg_exp) = names(short_names)
agg_exp=agg_exp[grep("_", names(agg_exp))]
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
  ) %>% 
  separate(region, into = c("cluster", "region"), sep = "_")

md_new = left_join(md_new, meta.data)
md_new = md_new %>% dplyr::filter(genotype!="WT")


res=pbmclapply(unique(md_new$cluster) %>% setNames(.,.), function(clst){  
  
  md = md_new %>% dplyr::filter(cluster==clst)
  counts1 = agg_exp2





counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]

maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
counts2 = counts1[maxs > 4, ]

exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
exp_gns = exp_gns[exp_gns > 0.8]

md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
md$region=factor(md$region)
md$sex=factor(md$sex)
md$age=factor(plyr::mapvalues(md$age, c("7 months", "14 months"), c("Y","O")))


  deseq = DESeqDataSetFromMatrix(
    countData = counts2[, c("gene", md$id)],
    colData = md,
    design =  ~ sex+age,
    tidy = T
  ) %>%
    estimateSizeFactors(type = "ratio")
  ############
  #sep
  ############
  dds_r = deseq[, deseq$region == "R"]
  dds_r$region=droplevels(dds_r$region)
  design(dds_r)=~t_s+t_c
  dds_r=dds_r %>%
    estimateDispersions(fitType="local") %>% 
    DESeq(test = "LRT",
          full = ~t_s+t_c,
          reduced = ~1,
          fitType = "local"
    )
  
  dds_c = deseq[, deseq$region == "C"]
  dds_c$region=droplevels(dds_c$region)
  design(dds_c)=~t_s+t_c
  dds_c=dds_c %>%
    estimateDispersions(fitType="local") %>% 
    DESeq(test = "LRT",
          full = ~ t_s+t_c,
          reduced = ~1,
          fitType = "local"
    )
  
  
  
  
  dds_m = deseq[, deseq$region == "M"]
  dds_m$region=droplevels(dds_m$region)
  design(dds_m)=~t_s+t_c
  dds_m=dds_m %>%
    estimateDispersions(fitType="local") %>% 
    DESeq(test = "LRT",
          full = ~ t_s+t_c,
          reduced = ~1,
          fitType = "local"
    )
  res_sep=left_join(left_join(results(dds_r, tidy=T) %>% dplyr::rename(gene=row, R=padj) %>% dplyr::select(gene,R),
               results(dds_m, tidy=T) %>% dplyr::rename(gene=row, M=padj) %>% dplyr::select(gene,M) ),
               results(dds_c, tidy=T) %>% dplyr::rename(gene=row, C=padj) %>% dplyr::select(gene,C)) 
  res_sep[is.na(res_sep)]=1
  
  #######
  #RM
  #######
  dds_rm = deseq[, deseq$region != "C"]
  dds_rm$region=droplevels(dds_rm$region)
  design(dds_rm)=~region+t_s+t_c
  dds_rm=dds_rm %>%
    estimateDispersions(fitType="local") %>% 
  DESeq(test = "LRT",
        full = ~ region+t_s+t_c,
      reduced = ~region,
      fitType = "local"
    )
  
  res_rm=results(dds_rm, tidy=T) %>% dplyr::rename(gene=row, RM=padj) %>%
    dplyr::select(gene, RM) %>% mutate(RM=replace_na(RM,1))
 res_rm=left_join(res_rm, res_sep)
 
 design(dds_rm)=~region+t_s+t_c+region:(t_s+t_c)
 dds_rm=dds_rm %>%
    estimateDispersions(fitType="local") %>% 
    DESeq(test = "LRT",
          full = ~ region+t_s+t_c+region:(t_s+t_c),
          reduced = ~region+t_s+t_c,
          fitType = "local"
    )
  
res_rm_int=left_join(left_join(results(dds_rm, tidy=T) %>% dplyr::rename(gene=row),
                    coefficients(dds_rm) %>%
                      as.data.frame() %>%
                      rownames_to_column("gene")), res_rm)

rm_fdr=res_rm_int %>% dplyr::filter(R<0.1|M<0.1| RM<0.1)
rm_fdr$fdr=p.adjust(rm_fdr$pvalue,method = "BH")
res_rm_int=left_join(res_rm_int,rm_fdr %>% dplyr::select(gene, fdr))%>%
  mutate(
    phi_M = atan2(t_s, t_c) %% (2 * pi),
    amp_M= sqrt(t_c ^ 2 + t_s ^ 2),
    phi_hr_M = (12 / pi) * phi_M, 
    phi_R = atan2((t_s+regionR.t_s), t_c+regionR.t_c) %% (2 * pi),
    amp_R= sqrt((t_c+regionR.t_c) ^ 2 + (t_s+regionR.t_s) ^ 2),
    phi_hr_R = (12 / pi) * phi_R, 
  ) 

#######
# RC
#######
# keep only Rostral and Caudal (drop Medial)
dds_rc = deseq[, deseq$region != "M"]
dds_rc$region = droplevels(dds_rc$region)

# main-effect test
design(dds_rc) = ~ region + t_s + t_c
dds_rc = dds_rc %>%
  estimateDispersions(fitType = "local") %>%
  DESeq(test    = "LRT",
        full    = ~ region + t_s + t_c,
        reduced = ~ region,
        fitType = "local")

res_rc = results(dds_rc, tidy = TRUE) %>%
  dplyr::rename(gene = row, RC = padj) %>%
  dplyr::select(gene, RC) %>%
  mutate(RC = replace_na(RC, 1))
res_rc = left_join(res_rc, res_sep)

# interaction test
design(dds_rc) = ~ region + t_s + t_c + region:(t_s + t_c)
dds_rc = dds_rc %>%
  estimateDispersions(fitType = "local") %>%
  DESeq(test    = "LRT",
        full    = ~ region + t_s + t_c + region:(t_s + t_c),
        reduced = ~ region + t_s + t_c,
        fitType = "local")

res_rc_int = left_join(
  left_join(
    results(dds_rc, tidy = TRUE) %>% dplyr::rename(gene = row),
    coefficients(dds_rc) %>%
      as.data.frame() %>%
      rownames_to_column("gene")
  ),
  res_rc
)

rc_fdr = res_rc_int %>%
  dplyr::filter(R < 0.1 | C < 0.1 | RC < 0.1) %>%
  mutate(fdr = p.adjust(pvalue, method = "BH"))

res_rc_int = left_join(res_rc_int, rc_fdr %>% dplyr::select(gene, fdr)) %>%
  mutate(
    phi_C = atan2(t_s, t_c) %% (2 * pi),
    amp_C= sqrt(t_c ^ 2 + t_s ^ 2),
    phi_hr_C = (12 / pi) * phi_C, 
    phi_R = atan2((t_s+regionR.t_s), t_c+regionR.t_c) %% (2 * pi),
    amp_R= sqrt((t_c+regionR.t_c) ^ 2 + (t_s+regionR.t_s) ^ 2),
    phi_hr_R = (12 / pi) * phi_R, 
  ) 


#######
# MC
#######
# keep only Medial and Caudal (drop Rostral)
dds_mc = deseq[, deseq$region != "R"]
dds_mc$region = droplevels(dds_mc$region)

# main-effect test
design(dds_mc) = ~ region + t_s + t_c

dds_mc = dds_mc %>%
  estimateDispersions(fitType = "local") %>%
  DESeq(test    = "LRT",
        full    = ~ region + t_s + t_c,
        reduced = ~ region,
        fitType = "local")

res_mc = results(dds_mc, tidy = TRUE) %>%
  dplyr::rename(gene = row, MC = padj) %>%
  dplyr::select(gene, MC) %>%
  mutate(MC = replace_na(MC, 1))
res_mc = left_join(res_mc, res_sep)

# interaction test
design(dds_mc) = ~ region + t_s + t_c + region:(t_s + t_c)
dds_mc = dds_mc %>%
  estimateDispersions(fitType = "local") %>%
  DESeq(test    = "LRT",
        full    = ~ region + t_s + t_c + region:(t_s + t_c),
        reduced = ~ region + t_s + t_c,
        fitType = "local")

res_mc_int = left_join(
  left_join(
    results(dds_mc, tidy = TRUE) %>% dplyr::rename(gene = row),
    coefficients(dds_mc) %>%
      as.data.frame() %>%
      rownames_to_column("gene")
  ),
  res_mc
)

mc_fdr = res_mc_int %>%
  dplyr::filter(M < 0.1 | C < 0.1 | MC < 0.1) %>%
  mutate(fdr = p.adjust(pvalue, method = "BH"))

res_mc_int = left_join(res_mc_int, mc_fdr %>% dplyr::select(gene, fdr)) %>%
  mutate(
    phi_C = atan2(t_s, t_c) %% (2 * pi),
    amp_C= sqrt(t_c ^ 2 + t_s ^ 2),
    phi_hr_C = (12 / pi) * phi_C, 
    phi_M = atan2((t_s+regionM.t_s), t_c+regionM.t_c) %% (2 * pi),
    amp_M= sqrt((t_c+regionM.t_c) ^ 2 + (t_s+regionM.t_s) ^ 2),
    phi_hr_M = (12 / pi) * phi_M, 
  ) 


res_int_all=list(RM=res_rm_int, RC=res_rc_int, MC=res_mc_int)

}, ignore.interactive = T)



saveRDS(res, 
        "~/desp1/precast/precast_final_with_ros_caud/analysis2/rhythmicity_interaction/res_interaction_APP.rds",
        compress = T)
