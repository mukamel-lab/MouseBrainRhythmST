library(tidyverse)
res_sig <- readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/downsample/coefs_wt_ds.rds") %>%
  dplyr::filter(padj<0.1)
res_sig$sig=ifelse(res_sig$padj<0.05, "FDR<0.05", "FDR<0.1")
res_sig$cluster=gsub("Olfactory Tubercle", "Amygdala", res_sig$cluster)
sig_colors=c("steelblue1","dodgerblue3")
  names(sig_colors)=c("FDR<0.05", "FDR<0.1")
  snms=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
  res_sig$cluster=plyr::mapvalues(res_sig$cluster, snms, names(snms))
res_sig$cluster=factor(res_sig$cluster, levels= names(snms))
  bar_sig=ggplot(res_sig  , aes(x=cluster, fill=sig))+
  geom_bar()+theme_bw()+
    scale_fill_manual(name="Significance Level",
                      values = sig_colors)+
  xlab("Cluster")+ylab("Rhythmic DEG Counts")+
    ggtitle("Number of Rhythmic Genes Per Cluster")+
  theme(title = element_text(size=18), axis.title =
          element_text(size=16), 
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text= element_text(size=15),
        axis.text.x = element_text(angle=70,hjust=1))

  
bar_sig+arial_theme+theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5 ), panel.grid.major.x = element_blank(),
                          panel.grid.minor.x = element_blank(),panel.grid.minor.y =element_blank(), legend.position = "none")
  

dta=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/spat_plots_generic/spat_data_030B.rds")

  dta=left_join(dta %>% mutate(cluster=plyr::mapvalues(anot, snms, names(snms))),
                res_sig %>% na.omit() %>% group_by(cluster) %>% 
                  tally())
  
  plt1=ggplot(dta)+
    geom_point(aes(y=-1*imagerow,x=imagecol, col=n), size=.5)+
    theme_void()+
    scale_color_gradient(name = "", 
                         high = "#16304A", 
                         low = "#66AEE0", 
                         trans="log10", breaks=c(0,10,50,200)) +
    ggtitle("")+coord_fixed()+
    theme(plot.title=element_text(size = 12, hjust = 0.5),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(.9, "lines")
    ) 

plt1
  plt1/bar_sig

