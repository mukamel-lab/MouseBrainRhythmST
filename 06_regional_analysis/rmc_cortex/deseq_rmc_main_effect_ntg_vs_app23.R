library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(ggh4x)

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
md_new = md_new %>% dplyr::filter(region!="M")


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
    design =  ~ sex+age+region,
    tidy = T
  ) %>%
    estimateSizeFactors(type = "ratio")
 
  dds_wt = deseq[, deseq$genotype=="WT"] %>%
    DESeq(fitType = "local")
  
  dds_app = deseq[, deseq$genotype!="WT"] %>%
    DESeq(fitType = "local")
  
  res=left_join(results(dds_wt, tidy=T) %>% dplyr::rename(gene=row),
                results(dds_app, tidy=T) %>% dplyr::rename(gene=row),
                by=c("gene"),
                suffix = c("_NTG", "_APP23")
                )
  
  res    

  }, ignore.interactive = T)



saveRDS(res, 
        "~/desp1/precast/precast_final_with_ros_caud/analysis2/main_effect/res_region_app23_ntg_joint.rds",
        compress = T)


res=res %>% bind_rows(.id="cluster")

p = ggplot(res, aes(x =log2FoldChange_NTG, y =log2FoldChange_APP23)) +
  facet_grid(.~cluster)+
  geom_abline(slope=1, lty=2)+
  sm_statCorr(corr_method = "spearman", color = "black",
              fullrange = TRUE)+
  geom_point(aes(x =log2FoldChange_NTG, y =log2FoldChange_APP23))+
  theme(text=element_text(size=6, family="ArialMT"))+
  coord_fixed()+
  xlab("NTG L2FC (Rostal/Caudal)")+
  ylab("APP23 L2FC (Rostal/Caudal)")

library(ggrastr)

p = ggplot(res%>% dplyr::filter(padj_NTG<0.05|padj_APP23<0.05), aes(x = log2FoldChange_NTG, y = log2FoldChange_APP23)) +
  facet_grid(. ~ cluster) +
  geom_abline(slope = 1, linetype = 2) +
  sm_statCorr(
    corr_method = "spearman",
    color = "black",
    fullrange = TRUE,
    text_size = 6 / .pt
      
  ) +
  rasterise(
    geom_point(size=.2),
    dpi = 600
  ) +
  theme(
    text = element_text(size = 6, family = "ArialMT"),
    strip.text = element_text(size = 6, family = "ArialMT"),
    axis.title = element_text(size = 6, family = "ArialMT"),
    axis.text = element_text(size = 6, family = "ArialMT")
  ) +
  coord_fixed() +
  xlab("NTG L2FC (Rostral/Caudal)") +
  ylab("APP23 L2FC (Rostral/Caudal)")

pdf("~/desp1/precast/precast_final_with_ros_caud/analysis2/main_effect/l2fc_Ros_caud_APPvsWT.pdf", height = 2.16, width = 5.45)
p
dev.off()
