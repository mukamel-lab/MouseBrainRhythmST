library(tidyverse)
library(patchwork)
library(pbmcapply)
library(ggh4x)
library(scales)

setwd("~/desp1/precast/precast_final_with_RMC_DV/")

counts=readRDS("analysis/deseq_rhyth/hip_regional/ncounts_NTG.rds")

coefs=readRDS("analysis/deseq_rhyth/hip_regional/res_NTG.rds") %>%
  
  mutate(
    phi = atan2(t_s, t_c) %% (2 * pi),
    amp = sqrt(t_c ^ 2 + t_s ^ 2),
    rel_amp = amp / log2(baseMean + 1),
    phi_hr = (12 / pi) * phi
  ) 


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


std_dev = function(x) {
  mean = mean(x)
  sd = sd(x)
  list(y=mean, ymin=mean-sd,
       ymax=mean+sd)
}
interact_res = readRDS("analysis/deseq_rhyth_dv_interact/hip_regional_NTG_reg_int.rds")%>% 
  dplyr::select(gene,age, padj)%>%
  mutate(sig=ifelse(padj<0.1,comparison,"ns" ))

clock_g = c("Arntl",
            "Cry1" ,
            "Cry2",
            "Dbp" ,
            "Per1" ,
            "Per2",
            "Per3",
            "Nr1d1",
            "Nr1d2",
            "Bhlhe41")

make_plot=function(gene_to_plot) {
  to_plt=counts %>% dplyr::filter(gene==gene_to_plot) 
  
  to_plt=rbind(to_plt, to_plt %>% 
                 mutate(time=plyr::mapvalues(time, paste0("ZT", 6*0:3 ),  
                                             paste0("ZT", 6*4:7 ))))
  to_plt$time=as.numeric(gsub("ZT", "",to_plt$time))
  
  
  crv=coefs %>% dplyr::filter(gene==gene_to_plot)
  
  crv=lapply(1:nrow(crv), function(i){
    mesor=to_plt %>%
      dplyr::filter(region==crv$region[i],
                    age==crv$age[i])  %>% 
      dplyr::select(sample, sex) %>%
      distinct()
    
    mesor=sum(mesor$sex=="F")/nrow(mesor)
    mesor=crv$Intercept[i]+crv$sex_M_vs_F[i]*mesor
    crv_spec=data.frame(time=seq(0,42,length.out=1000)) %>% 
      mutate(time_pi=(time %% 24)*2*pi/24)  %>%
      mutate(pred_exp=mesor+
               crv$t_c[i]*cos(time_pi)+
               crv$t_s[i]*sin(time_pi))%>%
      mutate(pred_exp=log2(2^pred_exp + 1))
    
    crv_spec$age=crv$age[i]
    crv_spec$region=crv$region[i]
    crv_spec
  }) %>% bind_rows()
  
  
  
  crv = crv %>%
    mutate(
      age = factor(age, levels=c("7 months", "14 months")),
      region = recode(region,
                  V="Ventral",
                  D="Dorsal"
                  
                      
      ),
      region=factor(region,levels=c("Dorsal","Ventral")))
  
  
  to_plt= to_plt%>%
    mutate(
      age = factor(age, levels=c("7 months", "14 months")),
      region = recode(region,
                      V="Ventral",
                      D="Dorsal"
                      
                      
      ),
      region=factor(region,levels=c("Dorsal","Ventral")))
  
  
  reg_col=ggsci::pal_d3()(5)[4:5] %>% c("Dorsal","Ventral")
  
  signif_df=interact_res %>%
    dplyr::filter(gene==gene_to_plot, padj<0.1) %>%
    group_by(age) %>%
    mutate(padj2=min(padj)) %>%
    mutate(lab=case_when(padj2<.01 ~ "***", 
                         padj2<0.05 ~ "**",.default = "*"
    )) %>%
    dplyr::select(age, lab) %>% distinct() %>%
    mutate(
      age = factor(age, levels=c("7 months", "14 months")))
  
  
  plt=ggplot()+
    annotate("rect", xmin = 0, xmax = 12, ymin = -Inf, ymax = Inf, fill = "#F6F18F", alpha = 0.2) +
    annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf, fill = "#606161", alpha = 0.2)+
    annotate("rect", xmin = 24, xmax = 36, ymin = -Inf, ymax = Inf, fill = "#F6F18F", alpha = 0.2) +
    annotate("rect", xmin = 36, xmax = 42, ymin = -Inf, ymax = Inf, fill = "#606161", alpha = 0.2)+
    geom_point(data=to_plt,
               aes(x=time, y=l2expr, col=region),
               stat="summary",fun=mean, size=.5)+
    
    stat_summary(data =to_plt,aes(x=time, y=l2expr, col=region),
                 fun.data = "std_dev",
                 geom = "errorbar",width=0.02,linewidth = 0.21)+
    
    geom_smooth(data=crv, aes(x=time, y=pred_exp, col=region),
                inherit.aes = F, linewidth = 0.21)+
    facet_nested(age~.,switch = "y")+
    scale_color_manual(name="",
                       values=reg_col)+
    ggtitle(gene_to_plot)+
    theme_bw(base_size = 8, base_family = "ArialMT")+
    theme(axis.title = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.key  = element_rect(fill="white"),
          strip.background = element_rect(fill = "white",
                                          color = "black", 
                                          linetype = "blank"),
          strip.placement = "outside",
          #strip.text = element_text(size = 6,family = "ArialMT"),
          plot.title=element_text(size = 6, hjust=.5,family = "ArialMT", face = "italic"),
          #strip.text.y = element_text(angle = 0, size = 8, family = "ArialMT"))+
          strip.text = element_blank())+
    scale_x_continuous(   breaks = 0:3*12, 
                          labels = rep(c(0, 12), 2))+
    guides(col=guide_legend(title="",
                            override.aes = list(fill = NA)),
           linetype = guide_legend(override.aes = list(fill = NA)))+
    scale_y_continuous(breaks = breaks_pretty(n = 2))+
    geom_text(data=signif_df ,
              aes(x = 40,
                  y =min(to_plt$l2expr),  
                  label = lab),
              size =4, family="ArialMT", fontface="plain")
  
  plt
}


plts=lapply(  c("Dbp" ,
               "Nr1d1",
               "Nr1d2",
               "Arntl"), make_plot)




pdf("figures/extra/sin_plots_sig_clock_DV.pdf", height=3, width = 1.8)
(wrap_plots(plts, nrow=2)+plot_layout(guides = "collect"))& theme(legend.position = "bottom")

dev.off()  
