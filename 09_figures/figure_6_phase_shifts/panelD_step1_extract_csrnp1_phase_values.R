library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(scico)

cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
temp_col=cluster_color



x=vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/full/results.csv.gz")
phase_relamp <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/old_young_sep/phase_relamp.rds")

y=phase_relamp %>% dplyr::filter(gene =="Csrnp1")

z=y %>% dplyr::select(genotype, cluster,gene,rel_amp,age, phi, amp, phi_hr) %>%
  pivot_wider(id_cols =c( gene, age, cluster),
              names_from = genotype, values_from = c(rel_amp, amp, phi_hr, phi) )

z$amp_dif=2*(z$amp_APP-z$amp_WT)
z$rel_amp_dif=(z$rel_amp_APP-z$rel_amp_WT)
z$rel_amp_l2fc=log2(z$rel_amp_APP/z$rel_amp_WT)
z$phi_dif=(12/pi)*atan2(sin(z$phi_APP-z$phi_WT),cos(z$phi_APP-z$phi_WT))

saveRDS(z %>% dplyr::select(cluster, age,phi_dif), "~/desp1/precast/prec_c25q25g3000/svg_plotting/csrnp1_phase/csrnp1_vals_new.rds")
