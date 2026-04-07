library(tidyverse)
library(pbmcapply)
library(lme4)
library(glue)
library(patchwork)
library(ggh4x)
library(ggrepel)
library(scales)

setwd("~/desp1/precast/precast_final_with_ros_caud/")

meta.data=readRDS("precast_metadata.rds")

NE=readRDS("~/desp1/precast/prec_c25q25g3000/objects/md_ser_neuroE.rds")

cluster_colors=readRDS("cluster_color.rds")

short_names = readRDS("short_names.rds")

cluster_colors=cluster_colors%>% mutate(sn=plyr::mapvalues(cluster, short_names, names(short_names)))

NE=NE %>% dplyr::select(sample, cell,predicted_activity)

meta.data=left_join(meta.data, NE)

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)),
                               time=factor(time, levels=paste0("ZT", 0:3*6))) 

md =meta.data %>% group_by(sample, time, t_s,
                           t_c,
                           sex, age, genotype, anot) %>% 
  summarise(neuroE=mean(predicted_activity), nspots=n())

md=split(md, md$anot)

md=lapply(md, function(df) { split(df, paste0(df$genotype,"_", df$age))})
res=pbmclapply(md, function(clst){
  pbmclapply(clst, function(df) {
    f1=lmer(neuroE ~1+t_s+t_c+(1|sex),data=df)
    f2=lmer(neuroE ~1+(1|sex),data=df)
    x=anova(f2,f1)
    p.val=x$`Pr(>Chisq)`[2]
    
    data.frame(
      cluster = df$anot[1],
      age = df$age[1],
      genotype = df$genotype[1],
      pval = p.val,
      t_s = mean(coefficients(f1)$sex[,2]),
      t_c = mean(coefficients(f1)$sex[,3]),
      intcp = mean(coefficients(f1)$sex[,1])
      
    ) %>%
      mutate(
        phi = atan2(t_s,t_c) %% (2 * pi),
        amp = sqrt(t_c ^ 2 + t_s ^ 2),
        rel_amp = amp / intcp
      )
    
  })  %>% bind_rows()}) %>% bind_rows(.id="cluster") %>%
  {df=.;rownames(df)=NULL;df}


res$padj=p.adjust(res$pval, method = "BH")


fdr_cut=0.05

p1=ggplot(res %>% mutate(clst2=ifelse(padj<fdr_cut, cluster, NA),
                        sig=ifelse(padj<fdr_cut, "sig", "ns"))%>%arrange(sig) %>%
            mutate(genotype=relevel(factor(genotype),"WT")),
          aes(y = rel_amp, x =12*phi/pi, col = clst2, size=sig))+
  geom_point(alpha=.7)+
  scale_size_discrete(range=c(2,4))+
  scale_color_manual(name="",
                     values = setNames(cluster_colors$cluster_color, cluster_colors$cluster))+
  theme_bw(base_size = 16) +
  guides(size="none")+
  ylab("Rel. Amplitude") +
  xlab("Phase") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  coord_polar(start =0)+
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24))+
  facet_grid(age~genotype, as.table = F)+
  labs(subtitle="NeuroEstimator score diurnal rhythmicity",
       title =  glue("Neuronal Activation"),
       caption=glue("Colored points indicate rhythmicity FDR<{fdr_cut}"))+
  theme(
    plot.title = element_text(
      angle = 0,
      size = 16,
      face = "bold",
      vjust = 1
    ),
    plot.subtitle = element_text(
      angle = 0,
      size = 14,
      face = "plain",
      vjust = 1
    ),
    plot.caption = element_text(size = 14, hjust =1, vjust = 0)
  )



p1+theme(legend.position = "none")



p2=ggplot(res %>% dplyr::filter(padj<0.1,grepl(" - ", cluster)) %>%
            mutate(genotype=relevel(factor(genotype),"WT"),
                   cluster=plyr::mapvalues(cluster, short_names, names(short_names))),
          aes(y = rel_amp, x =12*phi/pi, col = cluster))+
  geom_point(size=1)+
  geom_text_repel(aes(label = cluster))+
  scale_color_manual(name="",
                     values = setNames(cluster_colors$cluster_color, 
                                       cluster_colors$sn))+
  theme_bw(base_size = 16) +
  guides(size="none")+
  ylab("Rel. Amplitude") +
  xlab("Phase") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  coord_polar(start =0)+
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24))+
  facet_grid(age~genotype, as.table = F)+
  labs(subtitle="NeuroEstimator score diurnal rhythmicity",
       title =  glue("Neuronal Activation"),
       caption=glue("Colored points indicate rhythmicity FDR<{fdr_cut}"))+
  theme(
    plot.title = element_text(
      angle = 0,
      size = 16,
      face = "bold",
      vjust = 1
    ),
    plot.subtitle = element_text(
      angle = 0,
      size = 14,
      face = "plain",
      vjust = 1
    ),
    plot.caption = element_text(size = 14, hjust =1, vjust = 0)
  )



p2+theme(legend.position = "none")


p3=ggplot(res  %>%
            mutate(genotype=relevel(factor(genotype),"WT"),
                   cluster=plyr::mapvalues(cluster, short_names, names(short_names))),
          aes(y = rel_amp, x =12*phi/pi, col = cluster))+
  geom_point(size=1)+
  geom_text_repel(aes(label = cluster))+
  scale_color_manual(name="",
                     values = setNames(cluster_colors$cluster_color, 
                                       cluster_colors$sn))+
  theme_bw(base_size = 16) +
  guides(size="none")+
  ylab("Rel. Amplitude") +
  xlab("Phase") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
  coord_polar(start =0)+
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24))+
  facet_grid(age~genotype, as.table = F)+
  labs(subtitle="NeuroEstimator score diurnal rhythmicity",
       title =  glue("Neuronal Activation"),
       caption=glue("Colored points indicate rhythmicity FDR<{fdr_cut}"))+
  theme(
    plot.title = element_text(
      angle = 0,
      size = 16,
      face = "bold",
      vjust = 1
    ),
    plot.subtitle = element_text(
      angle = 0,
      size = 14,
      face = "plain",
      vjust = 1
    ),
    plot.caption = element_text(size = 14, hjust =1, vjust = 0)
  )



p3+theme(legend.position = "none")

res2=res %>%
  mutate(genotype=relevel(factor(genotype),"WT"),
         region=plyr::mapvalues(cluster, cluster_colors$cluster, 
                                as.character(cluster_colors$region))) %>%
  dplyr::filter(!grepl("Layer 4",cluster )) %>% 
  group_by(region, genotype, age) %>%
  summarise(t_s=mean(t_s), t_c=mean(t_c),intcp=mean(intcp)) %>%
  mutate(
    phi = atan2(t_s,t_c) %% (2 * pi),
    amp = sqrt(t_c ^ 2 + t_s ^ 2),
    phi_hr=12*phi/pi,
    rel_amp=amp/intcp,
    region=factor(region, levels=levels(cluster_colors$region))
  )

reg_col=cluster_colors %>% dplyr::select(region_color, region) %>% distinct()
reg_col$region_color[1:3]=ggsci::pal_d3()(3)
reg_col=setNames(reg_col$region_color, as.character(reg_col$region))

res2= res2 %>% ungroup()%>%
  mutate(genotype=plyr::mapvalues(as.character(genotype),
                                               c("WT","APP23"), c("NTG", "APP23")), 
                      genotype=factor(genotype, levels=c("NTG", "APP23")))
p4=ggplot(res2,aes(col=region,fill = region))+
  geom_segment(aes(x = 12*phi/pi, xend = 12*phi/pi, y = 0, 
                   yend = rel_amp),
               arrow = arrow(length = unit(0.03, "npc") ,
                             type   = "closed",
                             angle  = 40,), linewidth = .6,
               key_glyph = draw_key_point)+
  theme_minimal(base_size = 6, base_family="ArialMT")+
  guides(size="none")+
  ylab("") +
  xlab("") + 
  scale_color_manual(name="",
                     values=reg_col)+
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.2)+
  coord_polar(start =0)+
 
  scale_x_continuous(breaks = 0:3*6, limits = c(0, 24))+
  scale_y_continuous(breaks = pretty_breaks(n = 3))+
  facet_grid(age~genotype, as.table = F, switch = "y")+
  labs(title =  glue("Neuronal activation (NeuroEstimator score)"),
       subtitle=glue("Regional mean resultant vectors of diurnal rhythmicity"),
       caption="r=relative amplitude (0 - 0.23)\n\u03B8=acrophase (ZT hrs)")+
  theme(
    plot.title = element_text(
      angle = 0,
      family = "ArialMT",
      size = 7,
      face = "bold",
      vjust = 1
    ),
    plot.subtitle = element_text(
      angle = 0,
      size = 6,
      face = "plain",
      family = "ArialMT",
      vjust = 1
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
    plot.margin = margin(t = 1, r = 1, b = 0.2, l = 1, unit = "lines")

  )



p4 
