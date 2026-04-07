library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)

###
# change to current run directory
prec_run="~/desp1/precast/prec_run"
setwd(prec_run)

data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

seuInt=readRDS("seuInt.rds")
cols_cluster = pals::glasbey(length(unique(seuInt$cluster)))


samps=list(1:16,17:32,33:48,49:64,65)

######################
# add annotation
######################
annotation=data.frame(cluster=1:length(unique(seuInt$cluster)),
                      anotation=c("Cortex Layer 5/6a",
                                  "CA1/Other HPF",
                                  "Cortex Layer 6",
                                  "Dentate Gyrus",
                                  "Fiber Tracts 1",
                                  "Pallidum",
                                  "Cortex Layer 2/3",
                                  "CA3",
                                  "Striatum",
                                  "Cortex Layer 4",
                                  "Subiculum",
                                  "Fiber Tracts 2",
                                  "Cortical Amygdala",
                                  "Amygdala",
                                  "Lateral Ventricle",
                                  "Meninges",
                                  "Piriform Cortex"))



seuInt$cell=paste0(str_split_fixed(rownames(seuInt@meta.data), "-",2)[,1], "-1")
seuInt$anot=plyr::mapvalues(seuInt$cluster, annotation$cluster, annotation$anotation)



saveRDS(seuInt, "seuInt.rds",compress=F)
write.csv(annotation, "annotation.csv")

dir.create("annotated_plots")
setwd("annotated_plots")


for(cl in unique(seuInt$cluster)){
  col_sp_clust=rep("gray80", length(unique(seuInt$cluster)))
  col_sp_clust[as.numeric(cl)]=cols_cluster[as.numeric(cl)] 
  
  
  
  pdf(paste0(gsub(" ", "_", gsub("\\/", "-",annotation[cl,2] )), ".pdf"))
  for(smp in samps) {
  
    p_clust = SpaPlot(seuInt, batch = smp, point_size = .5,
                      cols = col_sp_clust  , combine = F )
    p_clust=lapply(1:length(smp), function(pl){p_clust[[pl]]+coord_fixed()+
      theme(legend.position = "none")+ggtitle(metadata$sample[smp[pl]])})
    
    p=(wrap_plots(p_clust) &NoLegend())+
    plot_annotation(title = annotation[cl,2]) & 
    theme(plot.title = element_text(hjust=0.5, face = "bold"))
  
  plot(p)}
  dev.off()
  
}


