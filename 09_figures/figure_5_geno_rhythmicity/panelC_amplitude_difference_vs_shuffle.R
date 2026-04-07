library(tidyverse)
library(glue)
library(pbmcapply)
library(dplyr)
library(grid)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/shuffle_amp")

res_young=pbmclapply(list.files("shuffle_young", full=T), 
                     function(i) {x=readRDS(i)
                     if (nrow(x)==1) {return(NULL)} else {x}
                    } ) %>% bind_rows(.id = "seed")

res_young_sum=res_young %>% group_by(cluster)%>%
  summarise(mean=mean(mean_dif), sd=sd(mean_dif))

true_young=readRDS("true_young.rds") 
true_young=left_join(true_young, res_young_sum)


res_old=pbmclapply(list.files("shuffle_old", full=T), readRDS) %>% bind_rows(.id = "seed")

res_old_sum=res_old %>% group_by(cluster)%>%
  summarise(mean=mean(mean_dif), sd=sd(mean_dif))

true_old=readRDS("true_old.rds") 
true_old=left_join(true_old, res_old_sum)

all=rbind(res_young %>% mutate(age="7 months"),
          res_old%>% mutate(age="14 months")) %>% mutate(shuf="shuffled")

all=bind_rows(list(all, true_old %>% mutate(age="14 months",shuf="unshuffled"),
          true_young %>% mutate(age="7 months",shuf="unshuffled")))

cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
reg_col=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
all$cluster=plyr::mapvalues(all$cluster, short_names, names(short_names))
all$cluster=factor(all$cluster, levels= rev(names(short_names))) 
all$age=factor(all$age, levels=c("7 months", "14 months"))

p <- ggplot(all, aes(x = cluster, y = mean_dif, col = shuf)) +
  stat_summary(fun = mean, geom = "point") +
  stat_summary(
    fun.y = mean,
    fun.ymin = function(x) mean(x) - 2*sd(x),
    fun.ymax = function(x) mean(x) + 2*sd(x),
    geom = "errorbar"
  ) +
  scale_y_continuous(limits = c(-.15, .15), oob = scales::squish) +
  coord_flip(clip = "off") +
  theme_minimal(base_size =8, base_family = "ArialMT") +
  xlab("") +
  geom_hline(yintercept = .7, linetype = 2) +
  geom_hline(yintercept = -.6, linetype = 2) +
  
  theme(legend.title = element_blank(), 
        legend.position = "right",
        plot.margin = unit(c(1, 3, 1, 1), "cm"),
        panel.spacing = unit(1.5, "cm"),
        plot.title = element_text(hjust=0.5),
        axis.title.y = element_text(lineheight = 5),
        axis.text.y = element_text(size=6, family = "ArialMT", color =rev(reg_col)))+
  ggtitle("Empirical Relative Amplitude Shuffle Distributions") +
  scale_color_manual(values = c("shuffled" = "grey45", "unshuffled" = "red3"),labels=c("shuffled", "data")) +
  facet_wrap(~ age)
# 
# emp_p=rbind(true_old %>% mutate(age="14 months",shuf="unshuffled"),
#             true_young %>% mutate(age="7 months",shuf="unshuffled")) %>%
#   ungroup()%>%
#   mutate(emp_p=2*pnorm(mean_dif, mean, sd, lower.tail=F)) %>%
#   group_by(age) %>% mutate(exp_fdr=p.adjust(emp_p, method="BH"))%>%
#   mutate(emp_p_lab=paste0(signif(exp_fdr,2)),
#          sig=case_when(exp_fdr<0.01~ "***", 
#                        exp_fdr<0.05~ "**",
#                        exp_fdr<0.1~"*",.default = ""))%>%
#   mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
#         cluster=factor(cluster, levels= rev(names(short_names))),
#           age=factor(age, levels=c("7 months", "14 months")))
# 
# emp_p=rbind(true_old %>% mutate(age="14 months",shuf="unshuffled"),
#             true_young %>% mutate(age="7 months",shuf="unshuffled")) %>%
#   ungroup()%>%
#   mutate(emp_p=2*pnorm(mean_dif, mean, sd, lower.tail=F)) %>%
#   group_by(age) %>% mutate(exp_fdr=p.adjust(emp_p, method="BH"))%>%
#   mutate(emp_p_lab=paste0(signif(exp_fdr,2)),
#          sig=case_when(exp_fdr<0.01~ "***", 
#                        exp_fdr<0.05~ "**",
#                        exp_fdr<0.1~"*",.default = ""))%>%
#   mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
#          cluster=factor(cluster, levels= rev(names(short_names))),
#          age=factor(age, levels=c("7 months", "14 months")))

emp_p=left_join(bind_rows(
  res_young %>% mutate(age = "7 months"),
  res_old   %>% mutate(age = "14 months")),
  
  
  rbind(true_old %>% mutate(age="14 months"),
            true_young %>% mutate(age="7 months")) %>% 
    dplyr::rename(tmd=mean_dif),by=c("cluster", "age")) %>%
  group_by(cluster, age)%>%
  summarise(emp_p=(sum(abs(mean_dif) >= abs(tmd)) + 1) / (n() + 1),
            .groups = "drop"
  ) %>%

  mutate(exp_fdr = p.adjust(emp_p, method = "BH")) %>%
  mutate(emp_p_lab=paste0(signif(exp_fdr,2)),
         sig=case_when(exp_fdr<0.01~ "***", 
                       exp_fdr<0.05~ "**",
                       exp_fdr<0.1~"*",.default = ""))%>%
  mutate(cluster=plyr::mapvalues(cluster, short_names, names(short_names)),
         cluster=factor(cluster, levels= rev(names(short_names))),
         age=factor(age, levels=c("7 months", "14 months")),
         shuf="",
         mean_dif=0)

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
    ymin = 0.16,
    ymax = 0.16,
    data=emp_p[i,]
  )
}
print(p)




