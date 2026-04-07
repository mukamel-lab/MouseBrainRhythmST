library(tidyverse)
library(Matrix)
library(glmnet)
library(pbmcapply)

setwd("~/desp1/precast/prec_c25q25g3000/")

data_dir="/cndd2/agelber/hal/qc_aligned"


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24))) %>% dplyr::filter(age=="7 months")

agg_c=readRDS("objects/deseq_norm_l2cpm.rds")

agg_c=lapply(agg_c, function(df) {df[,intersect(colnames(df), c("gene", meta.data$sample))]})

meta.data=split(meta.data, meta.data$genotype)


##############################################
################################################
res_sig = readRDS("deseq_rhyth2/WT_harmonic//par_im/res_sig.rds") %>%
  arrange(padj)

# res_sig = readRDS("deseq_rhyth2/both_no_geno_term/par_im/res_sig.rds") %>% dplyr::filter(padj<0.01) %>%
#   arrange(padj)

res_sig=split(res_sig, res_sig$cluster)


agg_c_wt=lapply(agg_c , function(df) {df[, intersect(colnames(df), meta.data$WT$sample)]})

agg_c_app=lapply(agg_c , function(df) {df[, intersect(colnames(df), meta.data$APP23$sample)]})

#####################################

predict_wt=pbmclapply(names(agg_c_wt), function(clst) {
  
ctx=agg_c_wt[[clst]] %>% t()


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


names(predict_wt)=names(agg_c_wt)

###########################################

predict_app=pbmclapply(names(agg_c_app), function(clst) {
  
  ctx=agg_c_app[[clst]] %>% t()
  
  
  ######################################################################
  
  ######################################################################
  
  
  ctx2_2_2=as.matrix(ctx[,res_sig[[clst]]$gene])
  
  md1=as.data.frame(meta.data$APP23)
  rownames(md1)=md1$sample
  
  md1=md1[rownames(ctx2_2_2),]
  
  
  z_c=predict_wt[[clst]][["z_c"]]
  
  
 
  
  z_s=predict_wt[[clst]][["z_s"]]
  
  
  
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

ggplot(both %>% dplyr::filter(cluster %in% levels(cluster)[1:11]),
       aes(x=cluster,y=circ_dif_hr,col=genotype))+
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(width = 0.9)) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "errorbar",
               position = position_dodge(width = 0.9))+
  coord_flip()+ggpubr::theme_pubclean()

