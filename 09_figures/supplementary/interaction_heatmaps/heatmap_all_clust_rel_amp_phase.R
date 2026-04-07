library(tidyverse)
library(ComplexHeatmap)
library(patchwork)



x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/full/results.csv.gz") %>%
  dplyr::filter(padj<0.05)
#  x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/young/results.csv.gz") 
phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/old_young_sep/phase_relamp.rds")

y=phase_relamp %>% dplyr::filter(gene %in% x$gene)

z=y %>% dplyr::select(genotype, cluster,gene,rel_amp,age, phi, amp, phi_hr) %>%
  pivot_wider(id_cols =c( gene, age, cluster),
              names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$rel_amp_dif=(z$rel_amp_APP-z$rel_amp_WT)
z$rel_amp_l2fc=log2(z$rel_amp_APP/z$rel_amp_WT)
z$phi_dif=(12/pi)*atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))
cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
temp_col=cluster_color

to_plt_ral2fc=z %>% 
  dplyr::select(gene, cluster, age,rel_amp_l2fc ) %>% 
  pivot_wider(id_cols = gene, names_from = c(cluster, age), values_from =rel_amp_l2fc )

clord=paste0(rep(names(cluster_color), each=2), rep(c("_young", "_old"), times=23))
rword=x %>%dplyr::select(cluster, gene) %>% distinct()

to_plt_ral2fc=to_plt_ral2fc %>% column_to_rownames("gene")

to_plt_ral2fc=as.matrix(to_plt_ral2fc)



to_plt_ral2fc=to_plt_ral2fc[rword$gene,clord]
Heatmap(to_plt_ral2fc, cluster_columns = F,cluster_rows = F)
############################
to_plt_ral2fc=z %>% 
  dplyr::select(gene, cluster, age,phi_dif ) %>% 
  pivot_wider(id_cols = gene, names_from = c(cluster, age), values_from =phi_dif )

clord=paste0(rep(names(cluster_color), each=2), rep(c("_young", "_old"), times=23))

to_plt_ral2fc=to_plt_ral2fc %>% column_to_rownames("gene")

to_plt_ral2fc=as.matrix(to_plt_ral2fc)



to_plt_ral2fc=to_plt_ral2fc[rword$gene,clord]

Heatmap(to_plt_ral2fc, cluster_columns = F, cluster_rows = F)
