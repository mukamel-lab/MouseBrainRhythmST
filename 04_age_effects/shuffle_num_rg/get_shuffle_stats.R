library(tidyverse)
library(glue)
library(pbmcapply)
library(dplyr)
library(grid)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/age_int_new/")

res_yo_rats=pbmclapply(list.files("shuffle_mean_amp_dif/shuffles/", full=T), 
                     readRDS) %>% bind_rows(.id="seed") %>% mutate(shuf="shuffled")


true_YO=readRDS("shuffle_mean_amp_dif/true_OY.rds") %>%  
  pivot_wider(id_cols = c("cluster"), names_from = "age", names_prefix = "amp_",
              values_from = starts_with("amp")) %>%
  mutate(dif=amp_O-amp_Y)


all=bind_rows(list(res_yo_rats, true_YO %>% mutate(shuf="unshuffled")))

cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
reg_col=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
all$cluster=plyr::mapvalues(all$cluster, short_names, names(short_names))
all$cluster=factor(all$cluster, levels= rev(names(short_names)))


p <- ggplot(all , aes(x = cluster, y = dif, col = shuf)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(
    fun.y = mean,
    fun.ymin = function(x) mean(x) - sd(x),
    fun.ymax = function(x) mean(x) + sd(x),
    geom = "errorbar"
  ) +
  #scale_y_continuous(limits = c(-.15, .15), oob = scales::squish) +
  coord_flip(clip = "off") +
  theme_minimal(base_size =6, base_family = "ArialMT") +
  xlab("") +
  theme(text=element_text(size=6, family = "ArialMT", color ="black"),
        legend.title = element_blank(), 
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        plot.margin = unit(c(1, 3, 1, 1), "cm"),
        panel.spacing = unit(1.5, "cm"),
        plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(size=6, family = "ArialMT", color =rev(reg_col)))+
  ylab("Log2(Num. rhythmic genes in 14 mo./7 mo.)") +
  scale_color_manual(values = c("shuffled" = "grey45", "unshuffled" = "red3"),labels=c("shuffled", "data"))

p
true_YO_rat2=readRDS("shuffle_num_rg/true_OY.rds")  %>%  
  pivot_wider(id_cols = c("cluster"), names_from = "age", 
              values_from = starts_with("n")) %>%
  mutate(true_n1_rat=log2(n1_O/n1_Y),
         true_n05_rat=log2(n05_O/n05_Y), 
         true_n01_rat=log2(n01_Y/n01_Y))


emp_p=left_join(res_yo_rats,true_MF_rat2, by="cluster") %>% group_by(cluster) %>% 
  summarise(emp_p=(sum(abs(n05_rat) >= abs(true_n05_rat)) + 1) / (n() + 1),
            .groups = "drop"
  ) %>%

  mutate(exp_fdr = p.adjust(emp_p, method = "BH")) %>%
  mutate(emp_p_lab=paste0(signif(exp_fdr,2)),
         sig=case_when(exp_fdr<0.01~ "***", 
                       exp_fdr<0.05~ "**",
                       exp_fdr<0.1~"*",.default = "")) %>%
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
         cluster=factor(cluster, levels= rev(names(short_names))),
         shuf="",
         n01_rat=0,
         n05_rat=0,
         n1_rat=0)

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
layer(data = data, stat = StatIdentity, position = PositionIdentity,
geom = ggplot2:::GeomCustomAnn,
inherit.aes = TRUE, params = list(grob = grob,
xmin = xmin, xmax = xmax,
ymin = ymin, ymax = ymax))
}


for(i in 1:nrow(emp_p)) {
  p <- p + annotation_custom2(
    grob = textGrob(emp_p$sig[i]
                    ),
    xmin = emp_p$cluster[i],
    xmax = emp_p$cluster[i],
    ymin = 3,
    ymax = 3,
    data=emp_p[i,]
  )
}
print(p+geom_hline(yintercept = 0))




