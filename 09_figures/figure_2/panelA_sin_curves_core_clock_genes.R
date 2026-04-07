library(tidyverse)
library(patchwork)
setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth")
cluster_order <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_order.rds")
cluster_color <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
coefs_counts_wt=readRDS("WT_harmonic/par_im/coefs_counts.rds")

#######################
########################


make_plt_df=function(gn, cl){
  ################################################
  to_pl_wt=coefs_counts_wt[[cl]] %>% dplyr::filter(gene==gn)
  to_pl_wt=rbind(to_pl_wt, to_pl_wt %>% 
                   mutate(time=plyr::mapvalues(time, paste0("ZT", 6*0:3 ),  
                                               paste0("ZT", 6*4:7 ))))
  to_pl_wt$time=as.numeric(gsub("ZT", "",to_pl_wt$time))
  to_pl_wt$cluster=cl
  
  curve_pl_wt=data.frame(time=seq(0,42,length.out=1000)) %>% mutate(time_pi=(time %% 24)*2*pi/24) %>%
    mutate(pred_exp=to_pl_wt$Intercept[1]+0.5*to_pl_wt$age_7.months_vs_14.months[1]+0.5*to_pl_wt$sex_M_vs_F[1]+to_pl_wt$a[1]*cos(time_pi)+to_pl_wt$b[1]*sin(time_pi))
  curve_pl_wt$cluster=cl
  
  return(list(dta=to_pl_wt, crv=curve_pl_wt))
  
}

plt=lapply(c("Rbm3"), function(gn){
dfs=lapply(cluster_order[c(1,11)], function(clst) {
  make_plt_df(gn,clst)
})
names(dfs)=cluster_order[c(1,11)]
dfs_dta=lapply(dfs, function(x) {
  x$dta
})
dfs_curves=lapply(dfs, function(x) {
  x$crv
})
dfs_dta=bind_rows(dfs_dta, .id="cluster")
dfs_curves=bind_rows(dfs_curves, .id="cluster")

p=ggplot()+geom_point(data=dfs_dta, aes(x=time, y=l2expr, col=cluster))+
  geom_smooth(data=dfs_curves, aes(x=time, y=pred_exp, col=cluster), inherit.aes = F)+
  theme_bw()+
  xlab("")+ylab("")+
  #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
  #ggtitle(paste0(gn, " in ", cl))++
  scale_color_manual(values =cluster_color )+
  scale_x_continuous(breaks=c(0, 6,12,18,24,30,36,42)) +
  guides(col="none")+
  # guides(col=guide_legend(title="Genotype",override.aes = list(fill = NA)),
  #        linetype = guide_legend(override.aes = list(fill = NA)))+
  theme(legend.key  = element_rect(fill="white"))+ggtitle(gn)

return(p)})




spec_col=cluster_color[c(1,11)]

plt_mean=lapply(cur$res %>% dplyr::filter(padj<0.1) %>% .[['row']], function(gn){
  dfs=lapply(cluster_order[c(1,11)], function(clst) {
    make_plt_df(gn,clst)
  })
  
  mean_logcpm = function(x) {
    log2(mean(2^x))
    
  }
  names(dfs)=cluster_order[c(1,11)]
  dfs_dta=lapply(dfs, function(x) {
    x$dta
  })
  dfs_curves=lapply(dfs, function(x) {
    x$crv
  })
  dfs_dta=bind_rows(dfs_dta, .id="cluster")
  dfs_curves=bind_rows(dfs_curves, .id="cluster")
  
  p=ggplot()+geom_point(data=dfs_dta, aes(x=time, y=l2expr, col=cluster),
                        stat="summary",fun=mean_logcpm)+
    
    stat_summary(data =dfs_dta,aes(x=time, y=l2expr, col=cluster),
                 fun.data = function(x) {
                   mean = mean_logcpm(x)
                   se = sd(x)/sqrt(length(x))
                   list(y=mean, ymin=mean-se,
                        ymax=mean+se)
                 },
                 geom = "errorbar",width=0.1)+
    geom_smooth(data=dfs_curves, aes(x=time, y=pred_exp, col=cluster), inherit.aes = F)+
    theme_bw()+
    xlab("")+ylab("")+
    #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
    #ggtitle(paste0(gn, " in ", cl))++
    scale_color_manual(values =spec_col )+
    scale_x_continuous(breaks=c(0, 6,12,18,24,30,36,42)) +
    #guides(col="none")+
    guides(col=guide_legend(title="Cluster",override.aes = list(fill = NA)),
           linetype = guide_legend(override.aes = list(fill = NA)))+
    theme(legend.key  = element_rect(fill="white"))+ggtitle(gn)
  
  
  return(p)})
wrap_plots(plt_mean)




plt_mean=lapply(cur$res %>% dplyr::filter(padj<0.1) %>% .[['row']], function(gn){
  dfs=lapply(cluster_order[c(1:5,11)], function(clst) {
    make_plt_df(gn,clst)
  })
  
  mean_logcpm = function(x) {
    log2(mean(2^x))
    
  }
  names(dfs)=cluster_order[c(1:5,11)]
  dfs_dta=lapply(dfs, function(x) {
    x$dta
  })
  dfs_curves=lapply(dfs, function(x) {
    x$crv
  })
  dfs_dta=bind_rows(dfs_dta, .id="cluster")
  dfs_curves=bind_rows(dfs_curves, .id="cluster")
  
  p=ggplot()+geom_point(data=dfs_dta, aes(x=time, y=l2expr, col=cluster),
                        stat="summary",fun=mean_logcpm)+
    
    stat_summary(data =dfs_dta,aes(x=time, y=l2expr, col=cluster),
                 fun.data = function(x) {
                   mean = mean_logcpm(x)
                   se = sd(x)/sqrt(length(x))
                   list(y=mean, ymin=mean-se,
                        ymax=mean+se)
                 },
                 geom = "errorbar",width=0.1)+
    geom_smooth(data=dfs_curves, aes(x=time, y=pred_exp, col=cluster), inherit.aes = F)+
    theme_bw()+
    xlab("")+ylab("")+
    #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
    #ggtitle(paste0(gn, " in ", cl))++
    scale_color_manual(values =cluster_color )+
    scale_x_continuous(breaks=c(0, 6,12,18,24,30,36,42)) +
    #guides(col="none")+
    guides(col=guide_legend(title="Cluster",override.aes = list(fill = NA)),
           linetype = guide_legend(override.aes = list(fill = NA)))+
    theme(legend.key  = element_rect(fill="white"))+ggtitle(gn)
  
  
  return(p)})
wrap_plots(plt_mean)




plt_mean=lapply(cur$res_sig$gene, function(gn){
  dfs=lapply(cluster_order[c(1,11)], function(clst) {
    make_plt_df(gn,clst)
  })
  
  mean_logcpm = function(x) {
    log2(mean(2^x))
    
  }
  names(dfs)=cluster_order[c(1,11)]
  dfs_dta=lapply(dfs, function(x) {
    x$dta
  })
  dfs_curves=lapply(dfs, function(x) {
    x$crv
  })
  dfs_dta=bind_rows(dfs_dta, .id="cluster")
  dfs_curves=bind_rows(dfs_curves, .id="cluster")
  
  p=ggplot()+geom_point(data=dfs_dta, aes(x=time, y=l2expr, col=cluster),
                        stat="summary",fun=mean_logcpm)+
    
    stat_summary(data =dfs_dta,aes(x=time, y=l2expr, col=cluster),
                 fun.data = function(x) {
                   mean = mean_logcpm(x)
                   se = sd(x)
                   list(y=mean, ymin=mean-se,
                        ymax=mean+se)
                 },
                 geom = "errorbar",width=0.1)+
    geom_smooth(data=dfs_curves, aes(x=time, y=pred_exp, col=cluster), inherit.aes = F)+
    theme_bw()+
    xlab("")+ylab("")+
    #xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
    #ggtitle(paste0(gn, " in ", cl))++
    scale_color_manual(values =spec_col )+
    scale_x_continuous(breaks=c(0, 6,12,18,24,30,36,42)) +
    #guides(col="none")+
    guides(col=guide_legend(title="Cluster",override.aes = list(fill = NA)),
           linetype = guide_legend(override.aes = list(fill = NA)))+
    theme(legend.key  = element_rect(fill="white"))+ggtitle(gn)
  
  
  return(p)})
wrap_plots(plt_mean[1:10])

