library(tidyverse)
library(glue)
library(ComplexHeatmap)
library(magrittr)


setwd("~/desp1/precast/prec_c25q25g3000/deseq_app/")

cluster_order = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")


res_sig = readRDS("res/res.rds") %>% filter(padj < 0.05) %>%
  mutate(cluster = gsub("Olfactory Tubercle", "Amygdala", cluster)) %>%
  dplyr::filter(cluster %in% cluster_order[1:11])
sig_tally = res_sig %>% group_by(gene) %>% tally()
genes_2 = sig_tally$gene[sig_tally$n >= 2]
res_sig = split(res_sig, res_sig$cluster) %>%
  lapply(., function(x) {
    intersect(x$gene, genes_2)
  })

names(res_sig)=gsub(" \\(L6b, CLA, EP\\)", "",
gsub("Dentate Gyrus",  "DG",
gsub("Cortex Layer ",
"CTX L", names(res_sig))))
cluster_order_shrt=gsub(" \\(L6b, CLA, EP\\)", "",
                        gsub("Dentate Gyrus",  "DG",
                             gsub("Cortex Layer ",
                                  "CTX L", cluster_order)))
res_sig_pl = make_comb_mat(res_sig)
res_sig_pl=res_sig_pl[comb_degree(res_sig_pl) >= 2]
res_sig_pl = res_sig_pl[comb_size(res_sig_pl) >= 4]
p_wt = UpSet(
  res_sig_pl,
  set_order = cluster_order_shrt[cluster_order_shrt %in% rownames(res_sig_pl)],
  comb_order = order(comb_degree(res_sig_pl), decreasing = T)
)
p_wt
