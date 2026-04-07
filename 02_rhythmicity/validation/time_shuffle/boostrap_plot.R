library(tidyverse)
setwd("/home/agelber/desp1/precast/prec_c25q25g3000/")
cluster_color = readRDS("objects/cluster_color.rds")

wt_res = readRDS("deseq_rhth/time_shuffle/bootstrap/wt_res.rds")
app_res = readRDS("deseq_rhth/time_shuffle/bootstrap/app_res.rds")

app_res2=data.frame(run=rep(paste0("s", 1:30), each=23),
                    cluster=rep(names(cluster_color), times=30)) 


app_res2=left_join(app_res2, app_res) %>% mutate(n=replace_na(n,0))


wt_res2=data.frame(run=rep(paste0("s", 1:30), each=23),
                    cluster=rep(names(cluster_color), times=30)) 


wt_res2=left_join(wt_res2, wt_res) %>% mutate(n=replace_na(n,0))


app_res2$genotype="APP23-TG"
wt_res2$genotype="NTG"
res=rbind(app_res2, wt_res2)
res2=res %>% group_by(genotype, cluster) %>% summarise(mean=mean(n), se=sd(n)/sqrt(n()))
genotype_colors=ggsci::pal_nejm()(2)
names(genotype_colors)=c("APP23-TG", "NTG")

ggplot(res2 %>% mutate(cluster=factor(cluster,
                                      levels=names(cluster_color))),
       aes(x=cluster, y=mean, fill=genotype,
                                 ymin = mean - se, ymax = mean + se))+
  geom_col(position=position_dodge(width = 0.6), width = 0.6)+
  geom_errorbar(position=position_dodge(width = 0.6), width=0.3 ) +
  theme_bw(base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))+
  scale_fill_manual(name="Genotype", values=genotype_colors)+xlab("")+ylab("Num. Rhythmic Genes")


