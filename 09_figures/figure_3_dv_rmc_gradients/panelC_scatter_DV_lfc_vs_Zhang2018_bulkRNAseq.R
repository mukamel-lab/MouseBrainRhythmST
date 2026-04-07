library(tidyverse)
library(ggsignif)
library(ggrepel)
setwd("~/desp1/precast/precast_final_with_RMC_DV/analysis2/hipp_seq")

arial_theme=readRDS("~/desp1/precast/precast_final_with_ros_caud/arial6_theme.rds")


base_theme = theme_bw(base_size = 6, base_family = "ArialMT") +
  arial_theme +
  theme(
    text=element_text(size = 6, family = "ArialMT", colour = "black"),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.position = "bottom", 
    plot.title = element_blank()
  )


hipposeq = vroom::vroom("elife-14997-supp3-v1.txt") %>% mutate(dv=ifelse(dg_v_fpkm>dg_d_fpkm, "Ventral","Dorsal"))

hip_DV_NTG <- readRDS("~/desp1/precast/precast_final_with_RMC_DV/analysis2/dv_main/hip_DV_NTG.rds") %>%
  dplyr::filter(cluster=="DGsg", row %in% hipposeq$gene_short_name) 


hip_DV_NTG$gene_short_name=hip_DV_NTG$row


hip_DV_NTG=left_join(hip_DV_NTG, hipposeq)
reg_col= c("Dorsal" = ggsci::pal_d3()(5)[4], "Ventral" = ggsci::pal_d3()(5)[5])

plt=ggplot(hip_DV_NTG, aes(x = dv, y = log2FoldChange)) +
  geom_point(aes(color = dv), size = 0.02, position = position_jitter(width = 0.2, height = 0)) +
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, outlier.shape = NA, width = 0.3) +
  labs(x = "Dorsal and Ventral markers in Cembrowski et. al.", y = "L2FC Dorsal vs. Ventral DGsg", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "right"
  )+
  geom_signif(
    comparisons = list(c("Dorsal", "Ventral")),
    map_signif_level = T,
    test = "t.test",
    test.args = list(paired = F),
    y_position = c(2), 
    tip_length = 0.02,
    size = 0.5,textsize = 6*5/14,family = "ArialMT"
  ) +
  scale_color_manual(name="", values = reg_col)+
  guides(color = guide_legend(override.aes = list(size = 1.5)))

# geom_text_repel(
#   data = to_lab,
#   aes(label = row),
#   force = 2,
#   max.overlaps = 50,
#   min.segment.length = .5
# )





pdf("~/desp1/precast/precast_final_with_RMC_DV/analysis2/hipp_seq/compare_dgsg_box_plot.pdf",
    height = 1.5,width = 2)
plt
dev.off() 

