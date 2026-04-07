library(tidyverse)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/downsample/spot10_res/")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

arial_theme=readRDS("~/desp1/precast/precast_final_with_ros_caud/arial6_theme.rds")


coefs_sig=  pbmclapply(list.files() %>% setNames(.,.),readRDS) %>% bind_rows(.id="run")

coefs_sig=coefs_sig %>% group_by(cluster, age, genotype) %>% summarise(mn=mean(n0.1), se=sd(n0.1)/sqrt(n())) 
coefs_sig=coefs_sig%>% ungroup()%>% mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
                              cluster=factor(cluster, levels=names(short_names)))
coefs_sig=coefs_sig%>% mutate(age=plyr::mapvalues(age, c("O","Y"), c("14 months", "7 months")))

p1= ggplot(coefs_sig %>% dplyr::filter(cluster %in% names(short_names)[1:11]), aes(x = cluster, y = mn, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mn - se, ymax = mn + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_grid2(age~., as.table = F, switch="y", axes = "all", remove_labels = "x")+
  scale_fill_manual(name = "Genotype", 
                    values = ggsci::pal_nejm()(2) %>% 
                      setNames(c("APP23","WT")),breaks=c("APP23","WT"),labels=c("APP23","NTG")
)+
  labs(x = NULL, y = "Number of rhythmic genes (FDR<0.1)") +
  theme_classic(base_size = 6, base_family = "ArialMT") +
  arial_theme +
  theme(
    text             = element_text(size = 6, family = "ArialMT", color = "black"),
    axis.text.y =  element_text(size = 6, family = "ArialMT", color = "black"),
    axis.text.x =  element_text(size = 6, family = "ArialMT", color = "black", angle = 90, vjust = 0.5, hjust = 1),
    legend.text      = element_text(size = 6, family = "ArialMT", color = "black"),
    legend.title     = element_text(size = 6, family = "ArialMT", color = "black"),
    legend.position  = "right",
    legend.key.size  = unit(0.4, "lines"),
    panel.grid       = element_blank(),
    strip.text = element_text(size = 6, family = "ArialMT", color = "black"),
    strip.placement = "outside",
    strip.background = element_blank(),
  )

p1
#################################################

setwd("~/desp1/precast/precast_final_with_ros_caud/analysis2/down_sample_num_rhyth/spots50_reg/res/")

coefs_sig=  pbmclapply(list.files() %>% setNames(.,.),readRDS) %>% bind_rows(.id="run")



coefs_sig= coefs_sig %>%
  group_by(age) %>%
  mutate(n=n/sum(n)) 

p2=ggplot(coefs_sig, aes(x = factor(age, levels=c("7 months", "14 months")),
                         y = n, fill=region)) +
  geom_col(width=0.4)+
  scale_fill_manual(
    name = "Region",
    values = ggsci::pal_d3()(3)
  ) +
  labs(
    x = "",
    y = "Proportion of rhythmic genes (FDR<0.1)"
  ) +
  theme_classic(base_size = 6, base_family = "ArialMT") +
  arial_theme +
  theme(
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
    legend.text = element_text(
      size = 6,
      family = "ArialMT",
      color = "black"
    ),
    panel.grid     = element_blank(),
    legend.position= "right",
    legend.key.size= unit(0.4, "lines"),
    
    legend.title   = element_text(
      size = 6,
      family = "ArialMT",
      color = "black"
    )
  )
pdf("~/desp1/precast/precast_final_with_ros_caud/figures2/bar_num_rhyth.pdf",height = 2, width =2.2)

(p1/p2)+plot_layout(guides = "collect")
dev.off()
