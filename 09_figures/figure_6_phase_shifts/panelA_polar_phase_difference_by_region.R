library(ggrastr)
library(tidyverse)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)


coefs=rbind(readRDS("young/interaction_results_readjusted_y.rds")  %>% 
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
         region=plyr::mapvalues(cluster, regions, names(regions)),
         region=factor(region, levels=names(regions)[!duplicated(names(regions))]),
         cluster=factor(cluster, levels=names(short_names)),
         rel_amp_dif=app_amp-wt_amp,
         phi_dif=atan2(sin(app_phi-wt_phi),cos(app_phi-wt_phi)) %% (2 * pi),
         phi_dif_hr=(12*phi_dif/pi+24)%%24,
         age="7 months"),
  readRDS("old/interaction_results_readjusted_o.rds")  %>% 
    mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
           region=plyr::mapvalues(cluster, regions, names(regions)),
           region=factor(region, levels=names(regions)[!duplicated(names(regions))]),
           cluster=factor(cluster, levels=names(short_names)),
           rel_amp_dif=app_amp-wt_amp,
           phi_dif=atan2(sin(app_phi-wt_phi),cos(app_phi-wt_phi)) %% (2 * pi),
           phi_dif_hr=(12*phi_dif/pi+24)%%24,
           age="14 months")
)%>% mutate(age=factor(age, levels=c("7 months", "14 months")))


df=coefs
amp_lims=max(abs(quantile(df$rel_amp_dif, .99)), quantile(df$rel_amp_dif, .99))

p=ggplot(df %>% dplyr::filter(fdr<0.1),
       aes(y = rel_amp_dif, x =phi_dif_hr)) +
  geom_point_rast(data=df  %>% dplyr::filter(fdr>=0.1,padj<.1),
                  aes(y = rel_amp_dif, x =phi_dif_hr),
                  col="grey70",alpha=0.2,
                  inherit.aes = F,size=.05)+

  geom_point(aes(col=cluster), size=.2) + 
  theme_minimal(base_size = 6, base_family="ArialMT")+
  guides(size="none")+
  ylab("") +
  xlab("") + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2, colour="black")+
  facet_grid2(age~region, switch = "y")+
  geom_hline(yintercept = 0,linetype = 2, linewidth = 0.2, colour = "black") + 
  coord_polar(start =0)+
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24))+
  scale_y_continuous(limits = c(-amp_lims,amp_lims), 
                     oob=scales::squish, 
                     breaks = pretty_breaks(n = 3))+
  ggsci::scale_color_igv(name="")+
  theme(
    plot.title = element_text(
      angle = 0,
      family = "ArialMT",
      size = 7,
      face = "bold",
      vjust = 1,
      hjust = .5
    ),
    
    legend.text = element_text(size = 6, family = "ArialMT"),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, family = "ArialMT"),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    plot.margin = margin(t = 1, r = 1, b = 0.2, l = 1, unit = "lines")
    
  )+guides(col=guide_legend(ncol = 1))

p

pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/polar_drg.pdf", height =6.38, width =7.44 )
plot(p)
dev.off()



