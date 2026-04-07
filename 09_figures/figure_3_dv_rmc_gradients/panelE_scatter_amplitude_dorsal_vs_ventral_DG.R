library(tidyverse)
library(ggrastr)
library(ggrepel)

setwd("~/desp1/precast/precast_final_with_RMC_DV/")

arial_theme = readRDS("arial6_theme.rds")

kegg_clock_genes_all = readRDS("~/desp1/precast/precast_final_with_ros_caud/kegg_clock_genes.rds")

counts = readRDS("analysis/deseq_rhyth/coefs_counts_NTG_joint_no_age_term.rds")

coefs=lapply(counts, function(x) {
  
  y=x$coefs
  colnames(y)=gsub("CA3so|CA3sp|DGmo|DGsg", "clst", colnames(y))
  y
  })
coefs=bind_rows(coefs, .id="cluster")

amps=coefs %>% dplyr::filter(Intercept>6 |gene %in% kegg_clock_genes_all) %>%mutate(
  V=sqrt((t_c+regionclst_V.t_c)^2+(t_s+regionclst_V.t_s)^2),
  D=sqrt((t_c)^2+(t_s)^2)
)

amps=amps %>% dplyr::select(gene,cluster,D,V)%>% 
                              pivot_longer(cols=c(-gene, -cluster), names_to = "group", values_to = "amp")
amps=amps%>% separate(group, into = c("region"), sep="_") %>% 
  pivot_wider(id_cols = c("gene","cluster"), names_from = "region", values_from = "amp")

amps=na.omit(amps)
short_names=readRDS("~/desp1/precast/precast_final_with_ros_caud/short_names.rds")

rhyth_res=readRDS("~/desp1/precast/precast_final_with_ros_caud/analysis/deseq_rhyth/WT/res_sig.rds")%>%
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)))
rhyth_res=rhyth_res %>% dplyr::select(gene, cluster,padj)

amps=left_join(amps,rhyth_res)
amps=amps %>% mutate(sig=case_when(padj<0.001~ "***",
                                   padj<0.01 ~ "**",
                                   padj<0.05 ~ "*", 
                                  .default = ""), padj=NULL)
genes_to_lab=c("Dbp","Per1", "Cirbp", "Hspa8", "Cry2", "Per3", "Rbm3", "Hspa5")
p1=ggplot(amps, aes(x=D, y=V)) +
  geom_abline(slope=1, intercept=0, color="darkred") +
  geom_point_rast(size=0.1,alpha=.2) +
  #geom_point_rast(data=filter(amps, sig=="***", !(gene %in% kegg_clock_genes_all)), color="red", size=0.3) +
  geom_point_rast(data=filter(amps, gene %in% kegg_clock_genes_all,D>.1& V>.1), color="blue", size=0.3) +
  geom_text_repel(
   # data=filter(amps, (gene %in% kegg_clock_genes_all & D>.2& V>.2)|sig=="***") %>% group_by(sig),
    data=filter(amps, gene %in% genes_to_lab),
    aes(label=gene),
    size=6*5/14,
    family="ArialMT",
    nudge_x = 0.05,         
    nudge_y = 0.05,
    color="black",
    force=2,
    max.overlaps=20,
    segment.size=0.2,
    min.segment.length=0
  ) +
  facet_grid(.~cluster, switch="y", as.table=F) +
  coord_fixed() +
  theme_bw(base_size=8, base_family="ArialMT") +
  arial_theme +
  theme(
    panel.grid.minor=element_blank(),
    legend.key=element_rect(fill="white"),
    strip.background=element_rect(fill="white", color="black", linetype="blank"),
    strip.placement="outside",
    strip.text=element_text(size=6, family="ArialMT"),
    axis.title=element_text(size=6, family="ArialMT"),
    strip.text.y=element_text(angle=0, size=6, family="ArialMT")
  ) +
  scale_x_continuous(limits=c(0,.75), oob=scales::squish) +
  scale_y_continuous(limits=c(0,.75), oob=scales::squish) +
  xlab("Dorsal") +
  ylab("Ventral")

pdf("figures/chapter2fig/fig2_2/scatter_rel_amp.pdf", height = 2.2,width = 4)
plot(
  p1
)
dev.off()
