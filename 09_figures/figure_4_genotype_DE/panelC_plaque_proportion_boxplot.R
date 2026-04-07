library(tidyverse)
setwd("/home/agelber/desp1/precast/prec_c25q25g3000")



md_ser_all=readRDS("figures/plaque_counts/plq_ano.rds")
  
sn=readRDS("objects/short_names2.rds")

md_ser_all$anot2=plyr::mapvalues(md_ser_all$anot, sn, names(sn))


md_samp=md_ser_all %>% group_by(sample, anot2, age, sex, genotype) %>%
  summarise(p_count=sum(plaque), s_count=n(), p_prop=p_count/s_count)

md_samp$anot2=factor(md_samp$anot2, levels=names(sn))
md_samp$age=factor(md_samp$age,levels=c("7 months","14 months"))

arial_theme=readRDS("~/desp1/precast/precast_final_with_ros_caud/arial6_theme.rds")

md_sum <- md_samp %>%
  group_by(anot2, genotype, age) %>%
  summarise(
    mu = mean(p_prop),
    se = sd(p_prop) / sqrt(n()),
    .groups = "drop"
  )

plt2 <- ggplot(md_sum,
               aes(
                 x = anot2,
                 y = mu,
                 ymin = mu - se,
                 ymax = mu + se
               )) +
  geom_col(
    position = position_dodge2(
      preserve = "single",
      padding  = 0.2,
      width    = 0.6
    ),
    width = 0.3,
    alpha=0.6
  ) +
  geom_errorbar(
    position = position_dodge2(
      preserve = "single",
      padding  = 0.2,
      width    = 0.6
    ),
    width = 0.3
  ) +
  geom_jitter(
    data = md_samp,
    inherit.aes = FALSE,
    aes(x = anot2, y = p_prop),
    position = position_jitter(width = 0.15, height = 0),
    size = 0.5,
    alpha = 1
  ) +
  facet_grid(age ~ genotype, scales = "free_x", drop = FALSE) +
  theme_bw(base_size = 6, base_family = "ArialMT") +
  arial_theme +
  theme(
    axis.ticks.x         = element_blank(),
    strip.background     = element_blank(),
    strip.text           = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    axis.text.x          = element_text(angle = 90, hjust = 1),
    axis.title.x         = element_blank(),
    panel.border         = element_blank(),
    axis.line.x          = element_line(color = "black"),
    axis.line.y          = element_line(color = "black")
  ) +
  ylab("Proportion of ST spots\nwith a plaque")+
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  scale_y_continuous(
    limits = c( NA,0.1),
    oob    = scales::squish
  ) 


pdf("~/desp1/precast/prec_c25q25g3000/figures/plaque_counts/plq_prop.pdf", height = 3,width = 4)
print(plt2)
dev.off()

