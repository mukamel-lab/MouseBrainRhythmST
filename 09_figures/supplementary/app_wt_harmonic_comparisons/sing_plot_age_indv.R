library(tidyverse)
library(patchwork)
setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth/old_young_sep")

coefs_counts_young=list(wt=readRDS("WT_young_harmonic/par_im/coefs_counts.rds"),
                        app=readRDS("APP_young_harmonic/par_im/coefs_counts.rds"))


coefs_counts_old=list(wt=readRDS("WT_old_harmonic/par_im/coefs_counts.rds"),
                      app=readRDS("APP_old_harmonic/par_im/coefs_counts.rds"))

coefs_counts_all=list(wt=readRDS("../WT_harmonic/par_im/coefs_counts.rds"),
                      app=readRDS("../APP_harmonic/par_im/coefs_counts.rds"))

mean_logcpm = function(x) {
  log2(mean(2^x))
  
}

se_logcpm = function(x) {
  mean = mean_logcpm(x)
  se = sd(x)/sqrt(length(x))
  list(y=mean, ymin=mean-se,
       ymax=mean+se)
}

std_err = function(x) {
  mean = mean(x)
  se = sd(x)/sqrt(length(x))
  list(y=mean, ymin=mean-se,
       ymax=mean+se)
}

make_plot=function(coefs_list, gn , cl, age){
  ################################################
  to_pl_wt=coefs_list[["wt"]][[cl]] %>% dplyr::filter(gene==gn)
  to_pl_wt=rbind(to_pl_wt, to_pl_wt %>% 
                   mutate(time=plyr::mapvalues(time, paste0("ZT", 6*0:3 ),  
                                               paste0("ZT", 6*4:7 ))))
  to_pl_wt$time=as.numeric(gsub("ZT", "",to_pl_wt$time))
  wt_mf_rat=sum(to_pl_wt$sex=="F")/nrow(to_pl_wt)
  wt_age_rat=sum(to_pl_wt$age=="14 months")/nrow(to_pl_wt)
  mesor=to_pl_wt %>% dplyr::select(sample,l2expr) %>% distinct()
  
  mesor_wt=mean(mesor$l2expr)
  
  curve_pl_wt=data.frame(time=seq(0,42,length.out=1000)) %>% mutate(time_pi=(time %% 24)*2*pi/24) %>%
    mutate(pred_exp=mesor_wt+
             to_pl_wt$a[1]*cos(time_pi)+
             to_pl_wt$b[1]*sin(time_pi))
  
  
  curve_pl_wt$age=age
  curve_pl_wt$geno="NTG"
  #########################################
  to_pl_app=coefs_list[["app"]][[cl]] %>% dplyr::filter(gene==gn)
  to_pl_app=rbind(to_pl_app, to_pl_app %>% 
                    mutate(time=plyr::mapvalues(time, paste0("ZT", 6*0:3 ),  
                                                paste0("ZT", 6*4:7 ))))
  to_pl_app$time=as.numeric(gsub("ZT", "",to_pl_app$time))
  
  app_mf_rat=sum(to_pl_app$sex=="F")/nrow(to_pl_app)
  app_age_rat=sum(to_pl_app$age=="7 months")/nrow(to_pl_app)
  mesor=to_pl_app %>% dplyr::select(sample, l2expr) %>% distinct()
  
  mesor_app=mean(mesor$l2expr)
  
  
  curve_pl_app=data.frame(time=seq(0,42,length.out=1000)) %>% mutate(time_pi=(time %% 24)*2*pi/24) %>%
    mutate(pred_exp=mesor_app+
             to_pl_app$a[1]*cos(time_pi)+
             to_pl_app$b[1]*sin(time_pi))
  
  curve_pl_app$age=age
  curve_pl_app$geno="APP23-TG"
  #############################################
  
  to_pl=rbind(to_pl_wt, to_pl_app) 
  curve_pl=rbind(curve_pl_wt, curve_pl_app)
  
  to_pl$geno=plyr::mapvalues(to_pl$genotype, c("WT", "APP23"), c("NTG", "APP23-TG"))
  
  ggplot()+geom_point(data=to_pl, aes(x=time, y=l2expr, col=geno),
                      stat="summary",fun=mean)+
    
    stat_summary(data =to_pl,aes(x=time, y=l2expr, col=geno),
                 fun.data = "std_err",
                 geom = "errorbar",width=0.1)+
    geom_smooth(data=curve_pl, aes(x=time, y=pred_exp, col=geno), inherit.aes = F)+
    theme_bw(base_size = 24)+xlab("Zeitgeber Time")+ylab("Log2 Normalized Expression")+
    # ggtitle(paste0(age))+
    ggtitle(NULL)+
    guides(col=guide_legend(title="Genotype",override.aes = list(fill = NA)),
           linetype = guide_legend(override.aes = list(fill = NA)))+
    theme(legend.key  = element_rect(fill="white"))+ggsci::scale_color_nejm()+
    scale_x_continuous(breaks = 0:7*6)
}



#################################################
#######################################################
setwd("../../figures/fig3_app_wt_harm/sin_plots/")

dir.create("old")
setwd("old")
res_sig=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_age_spec_g/res_old.rds") %>% dplyr::filter(padj<0.1)

temp_col<- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = names(temp_col)
names(short_names) = c(
  "CTX23",
  "CTX4",
  "CTX5a",
  "CTX5b",
  "CTX6a",
  "CTXsp",
  "CA1",
  "CA3so",
  "CA3sp",
  "DGmo",
  "DGsg",
  "RHP",
  "COAa",
  "COAp",
  "Piriform",
  "Caudoputamen",
  "STRv",
  "Amygdala",
  "GP",
  "RN",
  "LV",
  "FP",
  "Meninges"
)
res_sig$clust=plyr::mapvalues(res_sig$cluster, short_names, names(short_names))

plts=lapply(1:nrow(res_sig), function(i) {
  try({
    (make_plot(coefs_counts_old,
                  res_sig[i,"gene"], 
                  res_sig[i,"cluster"], 
                  "14 months") &theme(axis.title = element_blank()))+
      ggtitle(res_sig[i,"gene"])
  })
})


for(i in 1:nrow(res_sig)){
  if(class(plts[[i]])[1]!="try-error"){
    if(!dir.exists(res_sig[i,"clust"])){
      dir.create(res_sig[i,"clust"])
    }
    nm=paste0( res_sig[i,"clust"],
               "/",
               res_sig[i,"gene"],
               "_",
               gsub("\\.", "_",signif(res_sig[i,"padj"],3)),
               ".rds"
    )
    p1=plts[[i]]
    saveRDS(p1, nm, compress=F)
  }
  
}

#################################
setwd("..")
dir.create("young")
setwd("young")
res_sig=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_age_spec_g/res_young.rds") %>% dplyr::filter(padj<0.1)
res_sig$clust=plyr::mapvalues(res_sig$cluster, short_names, names(short_names))


plts=lapply(1:nrow(res_sig), function(i) {
  try({
    ((make_plot(coefs_counts_young,
                res_sig[i,"gene"], 
                res_sig[i,"cluster"], 
                "7 months")) &theme(axis.title = element_blank()))+
      ggtitle(res_sig[i,"gene"])
  })
})

for(i in 1:nrow(res_sig)){
  
    if(class(plts[[i]])[1]!="try-error"){
      if(!dir.exists(res_sig[i,"clust"])){
        dir.create(res_sig[i,"clust"])
      }
    nm=paste0( res_sig[i,"clust"],
               "/",
               res_sig[i,"gene"],
               "_",
               gsub("\\.", "_",signif(res_sig[i,"padj"],3)),
               ".rds"
    )
    p1=plts[[i]]
    saveRDS(p1, nm,compress=F)
  }
  
}
