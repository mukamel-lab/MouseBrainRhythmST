library(tidyverse)
library(ggrastr)
library(scattermore)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)



coefs=rbind(readRDS("young/interaction_results_readjusted_y.rds") %>% mutate(age="7 months"),
            readRDS("old/interaction_results_readjusted_o.rds") %>% mutate(age="14 months"))

res=coefs %>%  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
                      region=plyr::mapvalues(cluster, regions, names(regions)),
                      region=factor(region, levels=names(regions)[!duplicated(names(regions))]), 
                      age=factor(age,levels=c("7 months","14 months")))
summary_df = res %>% 
  group_by(region, age) %>% 
  summarise(
    mean_diff   = mean(app_amp - wt_amp),
    prop_o_gt_y = sum(app_amp > wt_amp)/n(),
  ) %>% 
  mutate(
    label = sprintf(
      "delt=%.2f\nP=%.2f",
      mean_diff,
      prop_o_gt_y
    ),
    
    x =.6,
    y = .05
  ) 


p=ggplot(res) +
  
  geom_scattermore(aes(x = wt_amp, y = app_amp, col = region),
                   pointsize   = 15.2,
                   alpha = .3,
                   pixels      = c(1500,1500),
                   interpolate =F) +
  geom_abline(slope = 1, intercept = 0, linewidth = .3,linetype=3 ) +
  facet_grid2( age~ region, strip = strip_vanilla(), switch="y") +
  scale_color_manual(values = region_color) +
  theme_classic(base_size = 6, base_family = "ArialMT") +
  scale_y_continuous(limits = c(0.05, 1.5),
                     oob = scales::squish,
                     trans = "log10") +
  scale_x_continuous(limits = c(0.05, 1.5),
                     oob = scales::squish,
                     trans = "log10") + coord_fixed() +
  geom_text(
    data = summary_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 6 * 5 / 14,
    family = "ArialMT",
    vjust = 0
  ) +
  xlab("NTG Log2 Amplitude") +
  ylab("APP23 Log2 Amplitude") +
  theme(
    legend.position = "none",
    text = element_text(
      size = 6,
      family = "ArialMT",
      color = "black"
    ),
    axis.text = element_text(
      size = 6,
      family = "ArialMT",
      color = "black"
    ),
    strip.text = element_text(
      size = 6,
      family = "ArialMT",
      color = "black"
    ),
    strip.background = element_blank(),
    strip.placement = "outside"
  )



pdf("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/figures/scatter_amp.pdf", height = 2.2, width =4.5 )
plot(p)
dev.off()

