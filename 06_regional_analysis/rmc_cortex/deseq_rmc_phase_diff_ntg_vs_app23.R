# Load libraries
library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(ggh4x)
library(ggrastr)

setwd("~/desp1/precast/precast_final_with_ros_caud/analysis2/rhythmicity_interaction/")

# Load results
wt = readRDS("res_interaction.rds")
app = readRDS("res_interaction_APP.rds")

# Function to extract phase difference
get_phase = function(x) {
  x$RC %>%
    dplyr::select(gene, fdr, phi_R, phi_C) %>%
    mutate(phase_diff = (12/pi) * atan2(sin(phi_R - phi_C), cos(phi_R - phi_C)))
}

# Combine NTG and APP23 phase differences
res = bind_rows(list(
  NTG = lapply(wt, get_phase) %>% bind_rows(.id = "cluster"),
  APP23 = lapply(app, get_phase) %>% bind_rows(.id = "cluster")
), .id = "genotype")

# Reshape to wide format for NTG vs APP23
res_wide = res %>%
  pivot_wider(id_cols = c(gene, cluster), names_from = genotype, values_from = c(fdr, phase_diff))

# Keep genes rhythmic in either group
res_w_r = res_wide %>% 
  dplyr::filter(fdr_NTG < 0.1 | fdr_APP23<0.1 )

# Label rhythmicity
res_w_r = res_w_r %>%
  mutate(Sig = case_when(
    fdr_NTG < 0.1 & fdr_APP23 < 0.1 ~ "Both",
    fdr_NTG < 0.1 ~ "NTG Only",
    fdr_APP23 < 0.1 ~ "APP23 Only"
  )) %>%
  mutate(Sig = factor(Sig, levels = rev(c("Both", "NTG Only", "APP23 Only")))) %>%
  arrange(Sig)

# Plot
p = ggplot(res_w_r , aes(x = phase_diff_NTG, y = phase_diff_APP23)) +
  facet_grid(. ~ Sig) +
  geom_abline(slope = 1, linetype = 2) +
  sm_statCorr(
    corr_method = "spearman",
    color = "black",
    fullrange = TRUE,
    text_size = 6 / .pt,
    show_text = T
  ) +
  rasterise(
    geom_point( size = 0.2),
    dpi = 600
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 6, family = "ArialMT"),
    strip.text = element_text(size = 6, family = "ArialMT"),
    axis.title = element_text(size = 6, family = "ArialMT"),
    axis.text = element_text(size = 6, family = "ArialMT")
  ) +
  coord_fixed() +
  xlab("NTG Phase Difference (Ros-Caud, hours)") +
  ylab("APP23 Phase Difference (Ros-Caud, hours)") +
  ggsci::scale_color_bmj()

# Save to PDF
pdf("~/desp1/precast/precast_final_with_ros_caud/analysis2/main_effect/phasediff_Ros_caud_APPvsWT.pdf", height = 2.16/1.2, width = 5.45/1.2)
p
dev.off()

r2=split(res_w_r,res_w_r$Sig)

circular::cor.circular((pi/12)*r2$`APP23 Only`$phase_diff_NTG,(pi/12)*r2$`APP23 Only`$phase_diff_APP23, test=T )
circular::cor.circular((pi/12)*r2$`Both`$phase_diff_NTG,(pi/12)*r2$`Both`$phase_diff_APP23, test=T )
circular::cor.circular((pi/12)*r2$`NTG Only`$phase_diff_NTG,(pi/12)*r2$`NTG Only`$phase_diff_APP23, test=T )
