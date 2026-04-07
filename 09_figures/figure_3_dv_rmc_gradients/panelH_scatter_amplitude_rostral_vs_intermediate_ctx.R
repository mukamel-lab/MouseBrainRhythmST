library(tidyverse)
library(ggrastr)
setwd("~/desp1/precast/precast_final_with_ros_caud/")

short_names = readRDS("~/desp1/precast/precast_final_with_ros_caud/short_names.rds")
arial_theme=readRDS("arial6_theme.rds")

short_namesdf = data.frame(full = short_names, short = names(short_names)) %>%
  separate(full, into = c("full", "region"), sep = " - ") %>%
  separate(short, into = c("short", "R"), sep = "_") %>%
  mutate(short = gsub("23", "2/3", short))

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

# coefs_y=readRDS("analysis2/rhythmicity_interaction/res_interaction_7mo.rds")
# coefs_o=readRDS("analysis2/rhythmicity_interaction/res_interaction_14mo.rds")

coefs=readRDS("analysis2/rhythmicity_interaction/res_interaction.rds")
clock_g=readRDS("kegg_clock_genes_all.rds")

# df1=rbind(coefs_y$L23$RM %>% mutate(age="7 months"),
#           coefs_o$L23$RM %>% mutate(age="14 months")) %>%
#   mutate(age=factor(age, levels=c("7 months", "14 months")))
df1=coefs$L23$RM

p1=ggplot(df1 %>% filter(!is.na(fdr))) +
  
  geom_point_rast(data = ~filter(.x, fdr >= 0.1 & !gene %in% clock_g),
             aes(x = amp_R, y = amp_M),
             color = "gray70", size = 0.02) +
  geom_abline(slope = 1, lty = 2, color = "darkred") +
  geom_point(data = ~filter(.x, fdr < 0.1 & !gene %in% clock_g),
             aes(x = amp_R, y = amp_M),
             color = "red", size = 0.2) +
  
  
  geom_point(data = ~filter(.x, gene %in% clock_g),
             aes(x = amp_R, y = amp_M),
             color = "blue", size = 0.2) +
  geom_text_repel(data = ~filter(.x, gene %in% clock_g),
                  aes(x =amp_R, y = amp_M, label = gene),
                  size = 6*5/14, family = "ArialMT", hjust = -0.2, vjust = 0.3, color = "black",
                  force=2,
                  max.overlaps = 10, 
                  segment.size = 0.07,
                  min.segment.length = .5) +
  coord_fixed() +
  labs(x = "Rostral", y = "Medial") +
  base_theme +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 6, family = "ArialMT", colour = "black"))+
  scale_x_continuous(limits =c(0, 1), oob = scales::squish, breaks = c(0,.5,1))+
  scale_y_continuous(limits =c(0, 1), oob = scales::squish, breaks = c(0,.5,1))





###########
# df2=rbind(coefs_y$L23$RC %>% mutate(age="7 months"),
#               coefs_o$L23$RC %>% mutate(age="14 months")) %>%
#   mutate(age=factor(age, levels=c("7 months", "14 months")))
df2=coefs$L23$RC
gtl1= df2 %>% dplyr::filter(fdr<0.1, !gene %in% clock_g ) %>% arrange(desc(abs(amp_R-amp_C))) %>% pull(gene) %>% .[2*(1:5)]

p2=ggplot(df2 %>% filter(!is.na(fdr))) +
  
  geom_point_rast(data = ~filter(.x, fdr >= 0.1 & !gene %in% clock_g),
             aes(x = amp_R, y = amp_C),
             color = "gray70", size = 0.02) +
  geom_abline(slope = 1, lty = 2, color = "darkred") +
  geom_point(data = ~filter(.x, fdr < 0.1 & !gene %in% clock_g),
             aes(x = amp_R, y = amp_C),
             color = "red", size = 0.2) +
  
  
  geom_point(data = ~filter(.x, gene %in% clock_g),
             aes(x = amp_R, y = amp_C),
             color = "blue", size = 0.2) +
  
  geom_text_repel(data = ~filter(.x, gene %in% c(clock_g,gtl1)),
                  aes(x =amp_R, y = amp_C, label = gene),
                  size = 6*5/14, family = "ArialMT", hjust = -0.2, vjust = 0.3, color = "black",
                  force=2,
                  nudge_x = 0.03,
                  nudge_y = 0.03,
                  max.overlaps = 10, 
                  segment.size = 0.07,
                  min.segment.length = .1) +
  coord_fixed() +
  labs(x = "Rostral", y = "Caudal") +
  base_theme +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 6, family = "ArialMT", colour = "black"))+
  scale_x_continuous(limits =c(0,1), oob = scales::squish, breaks = c(0,.5,1))+
  scale_y_continuous(limits =c(0, 1), oob = scales::squish, breaks = c(0,.5,1))

###########

# df3=rbind(coefs_y$L23$MC %>% mutate(age="7 months"),
#           coefs_o$L23$MC %>% mutate(age="14 months")) %>%
#   mutate(age=factor(age, levels=c("7 months", "14 months")))

df3=coefs$L23$MC

gtl2= df3 %>% dplyr::filter(fdr<0.1, !gene %in% clock_g ) %>% arrange(desc(abs(amp_C-amp_M))) %>% pull(gene) %>% .[2*(1:5)]

p3=ggplot(df3 %>% filter(!is.na(fdr))) +
  
  geom_point_rast(data = ~filter(.x, fdr >= 0.1 & !gene %in% clock_g),
             aes(x = amp_M, y = amp_C),
             color = "gray70", size = 0.02) +
  geom_abline(slope = 1, lty = 2, color = "darkred") +
  geom_point(data = ~filter(.x, fdr < 0.1 & !gene %in% clock_g),
             aes(x = amp_M, y = amp_C),
             color = "red", size = 0.2) +
  
  
  geom_point(data = ~filter(.x, gene %in% clock_g),
             aes(x = amp_M, y = amp_C),
             color = "blue", size = 0.2) +
  
  geom_text_repel(data = ~filter(.x, gene %in% c(clock_g,gtl2)),
                  aes(x =amp_M, y = amp_C, label = gene),
                  size = 6*5/14, family = "ArialMT", hjust = -0.2, vjust = 0.3, color = "black",
                  force=2,
                  nudge_x = 0.03,
                  nudge_y = 0.03,
                  max.overlaps = 10, 
                  segment.size = 0.07,
                  min.segment.length = .1 ) +
  coord_fixed() +
  labs(x = "Medial", y = "Caudal") +
  base_theme +
  theme(strip.placement = "outside",
        strip.text = element_text(size = 6, family = "ArialMT", colour = "black"))+
  scale_x_continuous(limits =c(0,1), oob = scales::squish, breaks = c(0,.5,1))+
  scale_y_continuous(limits =c(0, 1), oob = scales::squish, breaks = c(0,.5,1))



pdf("figures2/scatter_rel_amp_l23.pdf", height = 2.2,width =6.6)
plot(
  (p1 |  p3 | p2)&
    theme(
      plot.margin   = margin(0, 0, 0, 0),
      panel.spacing.y = unit(3, "pt"),
      panel.spacing.x = unit(3, "pt")
    ))
dev.off()
