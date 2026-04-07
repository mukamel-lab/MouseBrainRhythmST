library(pbmcapply)
library(MetaCycle)
library(tidyverse)

setwd("~/desp1/precast/prec_c25q25g3000/metacycle/bulk_cort_paper/")
wt=vroom::vroom("wt_bulk_cort.csv") 

meta2d(infile="wt_bulk_cort.csv", 
         outdir=".",
         filestyle="csv", 
       timepoints=as.numeric(gsub("\\.1", "", gsub("X", "",colnames(wt)[2:ncol(wt)])))
)

app=vroom::vroom("app_bulk_cort.csv") 

meta2d(infile="app_bulk_cort.csv", 
       outdir=".",
       filestyle="csv", 
       timepoints=as.numeric(gsub("\\.1", "", gsub("X", "",colnames(app)[2:ncol(app)])))
)


meta2d_wt_bulk_cort = vroom::vroom("meta2d_wt_bulk_cort.csv") 

meta2d_app_bulk_cort = vroom::vroom("meta2d_app_bulk_cort.csv") 

x=left_join(meta2d_wt_bulk_cort, meta2d_app_bulk_cort, by= c("CycID"), suffix = c(".wt", ".app"))
p2=ggplot(x)+geom_point(aes(x=meta2d_rAMP.wt, y=meta2d_rAMP.app))+geom_abline(slope = 1, intercept = 0)
p1=ggplot(x)+geom_point(aes(x=meta2d_AMP.wt, y=meta2d_AMP.app))+geom_abline(slope = 1, intercept = 0)+scale_y_log10()+scale_x_log10()

meta2d_wt_bulk_cort = vroom::vroom("meta2d_wt_bulk_cort.csv") %>% dplyr::filter(meta2d_BH.Q<0.05)
x=left_join(meta2d_wt_bulk_cort, meta2d_app_bulk_cort, by= c("CycID"), suffix = c(".wt", ".app"))
p3=ggplot(x)+geom_point(aes(x=meta2d_rAMP.wt, y=meta2d_rAMP.app))+geom_abline(slope = 1, intercept = 0)
p4=ggplot(x)+geom_point(aes(x=meta2d_AMP.wt, y=meta2d_AMP.app))+geom_abline(slope = 1, intercept = 0)+scale_y_log10()+scale_x_log10()

(p3+coord_fixed()|p4+coord_fixed())

(p2+coord_fixed()|p1+coord_fixed())



barplot(c(nrow(meta2d_wt_bulk_cort), nrow(meta2d_app_bulk_cort)), names.arg = c("WT", "APP"))
