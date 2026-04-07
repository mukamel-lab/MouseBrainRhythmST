library(tidyverse)
setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth")

coefs_counts_wt=readRDS("WT_harmonic/par_im/coefs_counts.rds")
coefs_counts_app=readRDS("APP_harmonic/par_im/coefs_counts.rds")

#######################
########################


make_plot=function(gn, cl){
################################################
to_pl_wt=coefs_counts_wt[[cl]] %>% dplyr::filter(gene==gn)
to_pl_wt=rbind(to_pl_wt, to_pl_wt %>% 
              mutate(time=plyr::mapvalues(time, paste0("ZT", 6*0:3 ),  
                                          paste0("ZT", 6*4:7 ))))
to_pl_wt$time=as.numeric(gsub("ZT", "",to_pl_wt$time))


curve_pl_wt=data.frame(time=seq(0,42,length.out=1000)) %>% mutate(time_pi=(time %% 24)*2*pi/24) %>%
  mutate(pred_exp=to_pl_wt$Intercept[1]+to_pl_wt$age_7.months_vs_14.months[1]+to_pl_wt$sex_M_vs_F[1]+to_pl_wt$a[1]*cos(time_pi)+to_pl_wt$b[1]*sin(time_pi))
curve_pl_wt$genotype="WT"
#########################################
to_pl_app=coefs_counts_app[[cl]] %>% dplyr::filter(gene==gn)
to_pl_app=rbind(to_pl_app, to_pl_app %>% 
                 mutate(time=plyr::mapvalues(time, paste0("ZT", 6*0:3 ),  
                                             paste0("ZT", 6*4:7 ))))
to_pl_app$time=as.numeric(gsub("ZT", "",to_pl_app$time))



curve_pl_app=data.frame(time=seq(0,42,length.out=1000)) %>% mutate(time_pi=(time %% 24)*2*pi/24) %>%
  mutate(pred_exp=to_pl_app$Intercept[1]+to_pl_app$age_7.months_vs_14.months[1]+to_pl_app$sex_M_vs_F[1]+to_pl_app$a[1]*cos(time_pi)+to_pl_app$b[1]*sin(time_pi))
curve_pl_app$genotype="APP23"
#############################################

to_pl=rbind(to_pl_wt, to_pl_app)
curve_pl=rbind(curve_pl_wt, curve_pl_app)

#return(list(to_pl=to_pl, curve_pl=curve_pl))
ggplot()+geom_point(data=to_pl, aes(x=time, y=l2expr, col=genotype))+
  geom_smooth(data=curve_pl, aes(x=time, y=pred_exp, col=genotype), inherit.aes = F)+
  theme_bw()+
  #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
  #ggtitle(paste0(gn, " in ", cl))+
  ggtitle(cl)+
  guides(col=guide_legend(title="Genotype",override.aes = list(fill = NA)),
linetype = guide_legend(override.aes = list(fill = NA)))+
  theme(legend.key  = element_rect(fill="white"))+ggsci::scale_color_nejm()+
  scale_x_continuous(breaks=c(0, 6,12,18,24,30,36,42)) +
  theme(axis.text = element_text(size=14))
}

p1=make_plot("Dusp6","CA1")
p1=make_plot("Ddn","Dentate Gyrus-sg")

dbps=lapply(cluster_order, function(x){
  make_plot("Ptma", x)
})
names(dbps)=cluster_order
wrap_plots(dbps[cluster_order[c(1,20)]], ncol=1)+plot_layout(guides="collect")

names(cluster_order)=cluster_order

ant1=lapply(cluster_order[7:12], function(x){
  make_plot("Atn1", x)
})


ant1_to_pl=lapply(ant1, function(x) {x$to_pl}) %>% bind_rows(.id="cluster")
ant1_curve=lapply(ant1, function(x) {x$curve_pl}) %>% bind_rows(.id="cluster")


ggplot()+geom_point(data=ant1_to_pl, aes(x=time, y=l2expr, col=cluster))+
  geom_smooth(data=ant1_curve, aes(x=time, y=pred_exp, col=cluster), inherit.aes = F)+
  theme_bw()+
  #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
  #ggtitle(paste0(gn, " in ", cl))+
  # ggtitle(cl)+
  guides(col=guide_legend(title="Genotype",override.aes = list(fill = NA)),
         linetype = guide_legend(override.aes = list(fill = NA)))+
  theme(legend.key  = element_rect(fill="white"))+ggsci::scale_color_nejm()+facet_wrap(~genotype)



ggplot()+geom_point(data=ant1_to_pl %>% filter(cluster %in% cluster_order[c(7,8,9,11)]), aes(x=time, y=l2expr, col=cluster))+
  geom_smooth(data=ant1_curve %>% filter(cluster %in% cluster_order[c(7,8,9,11)]), aes(x=time, y=pred_exp, col=cluster), inherit.aes = F)+
  theme_bw()+
  #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
  #ggtitle(paste0(gn, " in ", cl))+
  # ggtitle(cl)+
  guides(col=guide_legend(title="Cluster",override.aes = list(fill = NA)),
         linetype = guide_legend(override.aes = list(fill = NA)))+
  theme(legend.key  = element_rect(fill="white"))+ggsci::scale_color_nejm()+facet_wrap(~genotype)

ggplot()+geom_point(data=ant1_to_pl %>% filter(cluster %in% cluster_order[c(1:1)]), aes(x=time, y=l2expr, col=genotype))+
  geom_smooth(data=ant1_curve %>% filter(cluster %in% cluster_order[c(7,8,9,11)]), aes(x=time, y=pred_exp, col=genotype), inherit.aes = F)+
  theme_bw()+
  #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
  #ggtitle(paste0(gn, " in ", cl))+
  # ggtitle(cl)+
  guides(col=guide_legend(title="Cluster",override.aes = list(fill = NA)),
         linetype = guide_legend(override.aes = list(fill = NA)))+
  theme(legend.key  = element_rect(fill="white"))+ggsci::scale_color_nejm()+facet_wrap(~cluster)


ggplot()+geom_point(data=ant1_to_pl, aes(x=time, y=l2expr, col=genotype))+
  geom_smooth(data=ant1_curve , aes(x=time, y=pred_exp, col=genotype), inherit.aes = F)+
  theme_bw()+
  #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
  #ggtitle(paste0(gn, " in ", cl))+
  # ggtitle(cl)+
  guides(col=guide_legend(title="Cluster",override.aes = list(fill = NA)),
         linetype = guide_legend(override.aes = list(fill = NA)))+
  theme(legend.key  = element_rect(fill="white"))+ggsci::scale_color_nejm()+facet_wrap(~cluster)
