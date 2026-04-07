library(tidyverse)
library(ggsignif)
library(ggrepel)
setwd("~/desp1/precast/precast_final_with_ros_caud/")

arial_theme=readRDS("arial6_theme.rds")


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


tasic = vroom::vroom("analysis/deseq_regional_DE_ros_v_caud/NIHMS1001532-DE_results.csv") %>%
  dplyr::filter(grepl("^L([2-9])", `VISp Cluster`)|grepl("^L([2-9])", `ALM Cluster`))

ros=str_split(tasic$`DEGenes\nHigher in ALM`, " ") %>% unlist(recursive = T) %>% na.omit() %>% unique() 
caud=str_split(tasic$`DEGenes\nHigher in VISp`, " ") %>% unlist(recursive = T) %>% na.omit() %>% unique()
ros1=setdiff(ros,caud)
caud1=setdiff(caud,ros)


res=readRDS("analysis/deseq_regional_DE_ros_v_caud/deseq_NTG.rds")
res=results(res, tidy=T, name="region_R_vs_C")
df=data.frame(row=c(ros1,caud1), du_in=c(rep("ALM", length(ros1)), rep("VISp", length(caud1))))

df2=left_join(df, res %>% dplyr::select(row, log2FoldChange, padj, baseMean)) %>% dplyr::filter(padj<0.1)
to_lab=rbind(df2%>% top_n(4,log2FoldChange), df2 %>% top_n(4,-1*log2FoldChange))
reg_col=ggsci::pal_d3()(3)[c(1,3)] %>% setNames(c("ALM", "VISp"))
plt=ggplot(df2, aes(x = du_in, y = log2FoldChange)) +
  geom_point(aes(color = du_in), size = 0.02, position = position_jitter(width = 0.2, height = 0)) +
  geom_boxplot(color = "black", fill = NA, linewidth = 0.3, outlier.shape = NA, width = 0.3) +
  labs(x = "Upregulated region in Tasic et. al.", y = "L2FC Rostral vs. Caudal", title = "") +
  base_theme +
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = 6, family = "ArialMT", colour = "black"),
    legend.position = "right"
  )+
  geom_signif(
    comparisons = list(c("ALM", "VISp")),
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



  

pdf("~/desp1/precast/precast_final_with_ros_caud/figures/extra/tasic_compare_box_plot.pdf",
    height = 1.5,width = 2)
plt
dev.off() 

