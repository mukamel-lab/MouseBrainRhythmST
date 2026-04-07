library(ComplexHeatmap)
library(tidyverse)

cluster_order <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/cluster_order.rds")

#################
#WT
##################
res_sig = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res.rds") %>%
  filter(padj<0.05)

sig_tally= res_sig %>% group_by(gene) %>% tally()
genes_2=sig_tally$gene[sig_tally$n>=3]
res_sig=split(res_sig, res_sig$cluster) %>%
  lapply(., function(x) {intersect(x$gene,genes_2)})

res_sig_pl=make_comb_mat(res_sig)
res_sig_pl=res_sig_pl[comb_degree(res_sig_pl) >= 3]
res_sig_pl=res_sig_pl[comb_size(res_sig_pl) >= 2]

p_wt=UpSet(res_sig_pl, set_order = cluster_order, 
           comb_order = order(comb_size(res_sig_pl), decreasing = T))

#################
#APP
##################

res_sig = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/res_sig.rds")

sig_tally= res_sig %>% group_by(gene) %>% tally()
genes_2=sig_tally$gene[sig_tally$n>=3]
res_sig=split(res_sig, res_sig$cluster) %>%
  lapply(., function(x) {intersect(x$gene,genes_2)})

res_sig_pl=make_comb_mat(res_sig)
res_sig_pl=res_sig_pl[comb_degree(res_sig_pl) >= 3]
res_sig_pl=res_sig_pl[comb_size(res_sig_pl) >= 2]

p_app=UpSet(res_sig_pl, set_order = cluster_order)



#################
#Both
##################

res_sig =  readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/both_harmonic/par_im/res.rds") %>%
  filter(padj<0.05)
sig_tally= res_sig %>% group_by(gene) %>% tally()
genes_2=sig_tally$gene[sig_tally$n>=3]
res_sig=split(res_sig, res_sig$cluster) %>%
  lapply(., function(x) {intersect(x$gene,genes_2)})

res_sig_pl=make_comb_mat(res_sig)
res_sig_pl=res_sig_pl[comb_degree(res_sig_pl) >= 3]
res_sig_pl=res_sig_pl[comb_size(res_sig_pl) >2]

p_both=UpSet(res_sig_pl, set_order = cluster_order)
