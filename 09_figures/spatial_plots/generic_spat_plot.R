library(tidyverse)

setwd("~/desp1/precast/prec_c25q25g3000/figures/spat_plots_generic")

dta=readRDS("spat_data_030B.rds")

##########################
# categorical
##########################
ggplot(dta)+ geom_point(aes(y=-1*imagerow,x=imagecol, col=anot))+
  theme_void()+
  guides(col="none")+
  ggtitle(glue("Slice {s}"))+
  scale_color_manual(values=cols_cluster)+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size=18))+
  coord_fixed()


############################
# continuous
############################

ggplot(dta)+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=gene), size=0.5)+
  theme_void()+
  scale_color_gradient(name="", low="gray90", high=cols,
                       limits=c(0,
                                quantile(dta$gene, .99)),
                       oob=scales::squish)+
  ggtitle("title")+coord_fixed()+
  theme(plot.title=element_text(size = 16, hjust = 0.5),
        legend.title=element_text(size = 14),
        legend.text = element_text(size = 12))