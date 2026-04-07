library(tidyverse)
library(scattermore)
library(ggh4x)
library(ggrastr)
library(ggrepel)
library(scales)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_reg_int2")

res= readRDS("~/desp1/precast/prec_c25q25g3000/deseq_reg_int2/regional_DRG_NTG.rds")%>% 
  bind_rows(.id="comparison")

res$amp_dif=res$r1_amp-res$r2_amp

res$phi_dif=(12/pi)*atan2(sin(res$r1_phi-res$r2_phi),
                             cos(res$r1_phi-res$r2_phi))

#x=res %>%dplyr::filter(comparison=="L23_DGsg")
#x=res %>%dplyr::filter(comparison=="DGsg_RN")
#x=res %>%dplyr::filter(comparison=="DGsg_FT")
x=res %>%dplyr::filter(comparison=="DGsg_CP")

p1=ggplot() +
  geom_point_rast(data=x %>% na.omit() %>% dplyr::filter(fdr>=0.1),
                  aes(y = amp_dif, x =phi_dif%%24),
                  col="grey60",alpha=0.3,
                  inherit.aes = F,size=.1)+
  geom_point(data=x %>% dplyr::filter(fdr<0.1,!(rhythmic_in_r1&rhythmic_in_r2)),
             aes(y = amp_dif, x =phi_dif%%24),col="steelblue", size=.3, inherit.aes = F) +
  geom_point(data=x %>% dplyr::filter(fdr<0.1,rhythmic_in_r1,rhythmic_in_r2),
             aes(y = amp_dif, x =phi_dif%%24),col="red", size=.5, inherit.aes = F) + 
  geom_text_repel(data=x %>% dplyr::filter(rhythmic_in_r1,rhythmic_in_r2,fdr<0.1)%>%
                    dplyr::slice_min(order_by = fdr,n = 20)%>%
                    dplyr::slice_min(order_by = abs(amp_dif),n = 10),
             aes(y = amp_dif, x =phi_dif%%24, label=gene), inherit.aes = F, min.segment.length = 0.0001) + 
  theme_minimal(base_size = 16) +
  
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.2, colour="black") + 
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2, colour="black")+
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24)) +
  scale_y_continuous(limits=c(-0.5,.5), oob=scales::squish, breaks = pretty_breaks(n=3))+
  coord_polar(start =0)+
  theme_bw(base_size = 6, base_family = "ArialMT")+
  theme(
    plot.title = element_text(
      angle = 0,
      family = "ArialMT",
      size = 6,
      face = "bold",
      vjust = 1,
      hjust = .5
    ),
    
    plot.caption = element_text(size = 7, family = "ArialMT",
                                hjust =0,  margin = margin(t = -2, b = 0, unit = "pt")),
    
    legend.text = element_text(size = 6, family = "ArialMT"),
    strip.background = element_blank(),
    strip.text = element_text(size = 6, family = "ArialMT"),
    strip.placement = "outside",
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6, family = "ArialMT",color="black"),
    panel.grid.minor.y = element_blank(),
    legend.key.size = unit(0.5, "lines"),
    plot.margin = margin(t = 1, r = 1, b = 0.2, l = 1, unit = "lines"),
    panel.grid.major = element_line(colour = "gray85")
    
  )

p1
pdf("~/desp1/precast/precast_final_with_ros_caud/figures/figure2/plots/polar_volcano_RMC.pdf",height = 2, width =6 )
plot(p1)
dev.off()
