library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(ggh4x)
library(ggrastr)

setwd("~/desp1/precast/precast_final_with_ros_caud/analysis2/rhythmicity_interaction/")

wt=readRDS("res_interaction.rds")
app=readRDS("res_interaction_APP.rds")
get_amps= function(x) {
  
  x$RC %>% dplyr::select(gene, fdr, amp_R, amp_C) %>%
    mutate(amp_diff=amp_R-amp_C)
}
res=bind_rows(list(
  NTG=lapply(wt, get_amps) %>% bind_rows(.id="cluster"),
  APP23=lapply(app, get_amps) %>% bind_rows(.id="cluster") 
), .id="genotype")


res_wide=res %>%
  pivot_wider(id_cols = c(gene, cluster), names_from = genotype, values_from = c(fdr, amp_diff))

res_w_r= res_wide %>% dplyr::filter(fdr_NTG<0.1 | fdr_APP23<0.1) 

res_w_r= res_w_r %>% mutate( Sig= case_when(fdr_NTG<0.1 & fdr_APP23<0.1 ~ "Both", 
                                                     fdr_NTG<0.1 ~ "NTG Only",
                                                     fdr_APP23<0.1 ~ "APP23 Only"))%>%
  mutate(Sig=factor(Sig, levels=rev(c("Both", "NTG Only",  "APP23 Only")))) %>% arrange(Sig)

p = ggplot(res_w_r, aes(x = amp_diff_NTG, y = amp_diff_APP23)) +
  facet_grid(. ~ cluster) +
  geom_abline(slope = 1, linetype = 2) +
  sm_statCorr(
    corr_method = "spearman",
    color = "black",
    fullrange = TRUE,
    text_size = 6 / .pt
      
  ) +
  rasterise(
    geom_point(aes(x = amp_diff_NTG, y = amp_diff_APP23, col=Sig),inherit.aes = F, size=.2),
    dpi = 600
  ) +
  theme(legend.position = "bottom",
    text = element_text(size = 6, family = "ArialMT"),
    strip.text = element_text(size = 6, family = "ArialMT"),
    axis.title = element_text(size = 6, family = "ArialMT"),
    axis.text = element_text(size = 6, family = "ArialMT")
  ) +
  coord_fixed() +
  xlab("NTG Log2 Amplitude (Ros-Caud)") +
  ylab("APP23 Log2 Amplitude (Ros-Caud)")+ggsci::scale_color_bmj()

p

pdf("~/desp1/precast/precast_final_with_ros_caud/analysis2/main_effect/ampdiff_Ros_caud_APPvsWT.pdf", height = 2.16, width = 5.45)
p
dev.off()
