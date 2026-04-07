library(tidyverse)
library(Matrix)
library(glmnet)
library(pbmcapply)

setwd("~/desp1/precast/prec_c25q25g3000/")

data_dir="/cndd2/agelber/hal/qc_aligned"

cluster_order=readRDS("objects/cluster_order.rds")

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


meta.data=split(meta.data, meta.data$genotype)

agg_c=readRDS("objects/deseq_norm_l2cpm.rds")

##############################################
################################################
res_sig = readRDS("deseq_rhyth2/WT_harmonic/par_im/res_sig.rds") %>%
  arrange(padj)

# res_sig = readRDS("deseq_rhyth2/both_no_geno_term/par_im/res_sig.rds") %>% dplyr::filter(padj<0.01) %>%
#   arrange(padj)

res_sig=split(res_sig, res_sig$cluster)


# res_sig=lapply(res_sig, function(df) {df[1:20,]})

agg_c_wt=lapply(agg_c , function(df) {df[, intersect(colnames(df), meta.data$WT$sample)]})
agg_c_wt=lapply(agg_c_wt , function(df) {
  df_t=sample(colnames(df), round(0.4*ncol(df)))
  df_tr=df[,df_t]
  df_tst=df[,setdiff(colnames(df), df_t)]
  return(list(train=df_tr, test=df_tst))
  })


agg_c_wt_train=lapply(agg_c_wt , function(df) {df$train})
agg_c_wt_test=lapply(agg_c_wt , function(df) {df$test})

agg_c_app=lapply(agg_c , function(df) {df[, intersect(colnames(df), meta.data$APP23$sample)]})

#####################################

train_wt=pbmclapply(names(agg_c_wt_train), function(clst) {
  
ctx=agg_c_wt_train[[clst]] %>% t()


######################################################################

######################################################################


ctx2_2_2=as.matrix(ctx[,res_sig[[clst]]$gene])

md1=as.data.frame(meta.data$WT)
rownames(md1)=md1$sample

md1=md1[rownames(ctx2_2_2),]

z_c= cv.glmnet(ctx2_2_2, md1$t_c,intercept=F )



z_c=as.matrix(coef(z_c)[2:(ncol(ctx2_2_2)+1)])


z_s = cv.glmnet(ctx2_2_2, md1$t_s, intercept=F )

z_s=as.matrix(coef(z_s)[2:(ncol(ctx2_2_2)+1)])



fit_s=ctx2_2_2 %*% z_s
fit_c=ctx2_2_2 %*% z_c

fit_=cbind(fit_c, fit_s)
fit_=as.data.frame(fit_)
colnames(fit_)=c("cos_pred","sin_pred")

fit_ =fit_ %>% rownames_to_column("samp")
fit_$t_c=md1$t_c
fit_$t_s=md1$t_s
fit_$time=md1$time
fit_$tr_phi=atan2(fit_$t_s, fit_$t_c) %%(2*pi)
fit_$tr=12*(atan2(fit_$t_s, fit_$t_c) %%(2*pi))/pi
fit_$pred=12*(atan2(fit_$sin_pred, fit_$cos_pred) %%(2*pi))/pi

fit_$pred_phi=atan2(fit_$sin_pred, fit_$cos_pred) %% (2*pi)

fit_$circ_dif=atan2(sin(fit_$pred_phi-fit_$tr_phi),cos(fit_$pred_phi-fit_$tr_phi))

fit_$circ_dif_hr=12*fit_$circ_dif/pi


mse=sum(fit_$circ_dif_hr^2)/nrow(fit_)
return(list(z_s=z_s, z_c=z_c, fit_=fit_, mse=mse))

})


names(train_wt)=names(agg_c_wt)

###########################################

predict_wt=pbmclapply(names(agg_c_wt_test), function(clst) {
  
  ctx=agg_c_wt_test[[clst]] %>% t()
  
  
  ######################################################################
  
  ######################################################################
  
  
  ctx2_2_2=as.matrix(ctx[,res_sig[[clst]]$gene])
  
  md1=as.data.frame(meta.data$WT)
  rownames(md1)=md1$sample
  
  md1=md1[rownames(ctx2_2_2),]
  
  
  z_c=train_wt[[clst]][["z_c"]]
  
  
  
  
  z_s=train_wt[[clst]][["z_s"]]
  
  
  
  fit_s=ctx2_2_2 %*% z_s
  fit_c=ctx2_2_2 %*% z_c
  
  fit_=cbind(fit_c, fit_s)
  fit_=as.data.frame(fit_)
  colnames(fit_)=c("cos_pred","sin_pred")
  
  fit_ =fit_ %>% rownames_to_column("samp")
  fit_$t_c=md1$t_c
  fit_$t_s=md1$t_s
  fit_$time=md1$time
  fit_$tr_phi=atan2(fit_$t_s, fit_$t_c) %%(2*pi)
  fit_$tr=12*(atan2(fit_$t_s, fit_$t_c) %%(2*pi))/pi
  fit_$pred=12*(atan2(fit_$sin_pred, fit_$cos_pred) %%(2*pi))/pi
  
  fit_$pred_phi=atan2(fit_$sin_pred, fit_$cos_pred) %% (2*pi)
  
  fit_$circ_dif=atan2(sin(fit_$pred_phi-fit_$tr_phi),cos(fit_$pred_phi-fit_$tr_phi))
  
  fit_$circ_dif_hr=12*fit_$circ_dif/pi
  
  
  mse=sum(fit_$circ_dif_hr^2)/nrow(fit_)
  return(list(z_s=z_s, z_c=z_c, fit_=fit_, mse=mse))
  
})


names(predict_wt)=names(agg_c_wt)

##############################################
###########################################

predict_app=pbmclapply(names(agg_c_app), function(clst) {
  
  ctx=agg_c_app[[clst]] %>% t()
  
  
  ######################################################################
  
  ######################################################################
  
  
  ctx2_2_2=as.matrix(ctx[,res_sig[[clst]]$gene])
  
  md1=as.data.frame(meta.data$APP23)
  rownames(md1)=md1$sample
  
  md1=md1[rownames(ctx2_2_2),]
  
  
  z_c=train_wt[[clst]][["z_c"]]
  
  
 
  
  z_s=train_wt[[clst]][["z_s"]]
  
  
  
  fit_s=ctx2_2_2 %*% z_s
  fit_c=ctx2_2_2 %*% z_c
  
  fit_=cbind(fit_c, fit_s)
  fit_=as.data.frame(fit_)
  colnames(fit_)=c("cos_pred","sin_pred")
  
  fit_ =fit_ %>% rownames_to_column("samp")
  fit_$t_c=md1$t_c
  fit_$t_s=md1$t_s
  fit_$time=md1$time
  fit_$tr_phi=atan2(fit_$t_s, fit_$t_c) %%(2*pi)
  fit_$tr=12*(atan2(fit_$t_s, fit_$t_c) %%(2*pi))/pi
  fit_$pred=12*(atan2(fit_$sin_pred, fit_$cos_pred) %%(2*pi))/pi
  
  fit_$pred_phi=atan2(fit_$sin_pred, fit_$cos_pred) %% (2*pi)
  
  fit_$circ_dif=atan2(sin(fit_$pred_phi-fit_$tr_phi),cos(fit_$pred_phi-fit_$tr_phi))
  
  fit_$circ_dif_hr=12*fit_$circ_dif/pi
  
  
  mse=sum(fit_$circ_dif_hr^2)/nrow(fit_)
  return(list(z_s=z_s, z_c=z_c, fit_=fit_, mse=mse))
  
})


names(predict_app)=names(agg_c_app)

##############################################


pred_sumry=data.frame(wt_mean_dif=unlist( sapply(predict_wt, function(df) {mean(df$fit_$circ_dif_hr)})),
                      app_mean_dif=unlist(sapply(predict_app, function(df) {mean(df$fit_$circ_dif_hr)})),
                      wt_mse=unlist( sapply(predict_wt, function(df) {mean(df$mse)})),
                      app_mse=unlist( sapply(predict_app, function(df) {mean(df$mse)}))
                      
) %>% 
  rownames_to_column("cluster")


wts=lapply(predict_wt, function(df) {df$fit_}) %>% bind_rows(.id="cluster") %>% mutate(genotype="WT")
apps=lapply(predict_app , function(df) {df$fit_}) %>% bind_rows(.id="cluster") %>% mutate(genotype="APP23")
both=rbind(wts, apps)
both$cluster=factor(both$cluster, levels = names(readRDS("objects/cluster_color.rds")))

p_vals=both %>% group_by(cluster) %>% 
  summarise(fdr=p.adjust(wilcox.test(circ_dif_hr[genotype=="WT"], 
                                     circ_dif_hr[genotype=="APP23"])[["p.value"]], 
                         method="BH"))

ggplot(both %>% dplyr::filter(cluster %in% levels(cluster)[1:11]),
       aes(x=cluster,y=circ_dif_hr,col=genotype))+
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "errorbar",
               position = position_dodge(width = 0.9))+
  geom_text(data=p_vals %>% 
              dplyr::filter(cluster %in% levels(cluster)[1:11]),
                                          aes(x=cluster,y=.2,
                                label=paste0("fdr: ", signif(fdr,2))), inherit.aes = F)+
  coord_flip()+ggpubr::theme_pubclean()+ggsci::scale_color_nejm()+ylab("ZT predicted - true")+
  theme(legend.position = "bottom")+ggtitle("Train on WT, Predict APP23")

both=both %>% mutate(time=factor(time, levels=paste0("ZT", 6*0:3)))
p_vals_zt=both %>% ungroup() %>% group_by(cluster, time) %>% 
  summarise(fdr=p.adjust(wilcox.test(circ_dif_hr[genotype=="WT"], 
                                     circ_dif_hr[genotype=="APP23"])[["p.value"]], 
                         method="BH"))

ggplot(both %>% dplyr::filter(cluster %in% levels(cluster)[1:11]),
       aes(x=cluster,y=circ_dif_hr,col=genotype))+
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "errorbar",
               position = position_dodge(width = 0.9))+
  geom_text(data=p_vals_zt %>% 
              dplyr::filter(cluster %in% levels(cluster)[1:11], fdr<0.1),
            aes(x=cluster,y=1,
                label=case_when(fdr<0.01~ "***", fdr<0.05~"**", fdr<0.1~"*")), 
            inherit.aes = F, size=7)+
  coord_flip()+ggpubr::theme_pubclean(base_size = 16)+ggsci::scale_color_nejm()+
  ylab("ZT predicted - true")+
  theme(legend.position = "bottom")+
  ggtitle("Train on WT, Predict APP23")+facet_wrap(~time)



predict_wt=lapply(predict_wt, function(lst) {lst$fit_}) %>% 
  bind_rows(.id="cluster") %>% mutate(genotype="WT")
predict_app=lapply(predict_app, function(lst) {lst$fit_}) %>% 
  bind_rows(.id="cluster") %>% mutate(genotype="APP23")

joint_pred=bind_rows(predict_wt, predict_app)

joint_pred=joint_pred %>% mutate(cluster=factor(cluster, levels=cluster_order))

ggplot(joint_pred  %>% dplyr::filter(cluster %in%p_vals$cluster[p_vals$fdr<0.05]),
       aes(x = tr_phi*12/pi, y = pred_phi*12/pi, col = genotype)) +
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.9)) +
  stat_summary(
    fun.y = mean,
    fun.ymin = function(x)
      mean(x) - sd(x) / sqrt(length(x)),
    fun.ymax = function(x)
      mean(x) + sd(x) / sqrt(length(x)),
    geom = "errorbar",
    position = position_dodge(width = 0.9),
    width=0.5
  ) +
 ggpubr::theme_pubclean() + ggsci::scale_color_nejm() +
  ylab("ZT predicted") +xlab("True ZT")+
  theme(legend.position = "bottom") + ggtitle("Train on WT, Predict APP23")
  


ggplot(both %>%
         dplyr::filter(cluster %in% levels(cluster)[1:11]),
       aes(x=cluster,y=circ_dif_hr,col=genotype))+
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "errorbar",
               width=0.4,
               position = position_dodge(width = 0.9))+
  geom_text(data=p_vals_zt %>% 
              dplyr::filter( fdr<0.1) %>% dplyr::filter(cluster %in% levels(cluster)[1:11])
              ,
            aes(x=cluster,y=1,
                label=case_when(fdr<0.01~ "***", fdr<0.05~"**", fdr<0.1~"*")), 
            inherit.aes = F, size=7)+
  coord_flip()+ggpubr::theme_pubclean(base_size = 16)+ggsci::scale_color_nejm()+
  ylab("ZT predicted - true")+
  theme(legend.position = "bottom")+
  ggtitle("Train on WT, Predict APP23")+facet_wrap(~time, nrow = 1)+ylab("")

