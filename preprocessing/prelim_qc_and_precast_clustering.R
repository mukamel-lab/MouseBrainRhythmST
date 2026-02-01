library(Seurat)
library(tidyverse)
library(customfuncs)
library(pals)
library(PRECAST)
library(glue)
library(pbmcapply)
library(patchwork)

setwd("~/desp1/precast/")
data_dir="/cndd2/agelber/hal/qc_aligned"

k_val=25
q_val=25
n_genes=3000
min_features=500

prec_run=glue("prec_c{k_val}q{q_val}g{n_genes}_with_cov")
dir.create(prec_run)

system(glue("cp ~/desp1/precast/prec_pipe_add_cov.R {prec_run}/invocation.R"))
system(glue("cp ~/desp1/precast/old_script/add_annotation.R {prec_run}"))

metadata=data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>%
  dplyr::filter(sample %in% list.files(data_dir))

ser_objs=pbmclapply(metadata$sample, function(s) {
  x=Load10X_Spatial(glue("{data_dir}/{s}/outs")) %>% subset(nFeature_Spatial > min_features)
  gns_to_use=rownames(x@assays$Spatial@counts)
  gns_to_use=gns_to_use[!grepl("Thy1|humanAPP", gns_to_use)]
  
  y=CreateSeuratObject(counts=x@assays$Spatial@counts[gns_to_use,])
  y$row=x@images$slice1@coordinates$imagecol
  y$col=-1*x@images$slice1@coordinates$imagerow
  y$sample=s
  y$age=factor(metadata$age[metadata$sample==y$sample[1]])
  y$time=factor(metadata$time[metadata$sample==y$sample[1]])
  y$sex=factor(metadata$sex[metadata$sample==y$sample[1]])
  y$genotype=factor(metadata$genotype[metadata$sample==y$sample[1]])
  return(y)
})
names(ser_objs)=metadata$sample

prec=CreatePRECASTObject(
  ser_objs,
  numCores_sparkx=15,
  gene.number=n_genes,
  selectGenesMethod="SPARK-X"
)

saveRDS(prec, glue("/home/AD/agelber/desp1/precast/{prec_run}/prec1.rds"), compress=F)

prec=AddAdjList(prec, platform="Visium")

prec=AddParSetting(
  prec,
  Sigma_equal=FALSE,
  verbose=TRUE,
  coreNum=15,
  int.model="EEE"
)

prec=PRECAST(prec, K=k_val, q=q_val)

prec=SelectModel(prec)

saveRDS(prec, glue("/home/AD/agelber/desp1/precast/{prec_run}/prec2.rds"), compress=F)

seuInt=IntegrateSpaData(
  prec,
  species="Mouse",
  covariates_use=c("age", "sex", "time", "genotype")
)

seuInt$sample=plyr::mapvalues(as.integer(seuInt$batch), from=1:65, to=metadata$sample)

seuInt=RunUMAP(
  seuInt,
  reduction="PRECAST",
  dims=1:ncol(seuInt@reductions$PRECAST@cell.embeddings)
)

saveRDS(seuInt, glue("/home/AD/agelber/desp1/precast/{prec_run}/seuInt.rds"), compress=F)

setwd(prec_run)

cols_cluster=pals::glasbey(length(unique(seuInt$cluster)))
p=SpaPlot(seuInt, batch=NULL, point_size=1.4, cols=cols_cluster, combine=F)
p=lapply(1:65, function(pl){p[[pl]]+coord_fixed()+theme(legend.position="none")+ggtitle(metadata$sample[pl])})

dir.create("slice_plots")
setwd("slice_plots/")
names(p)=metadata$sample

for(i in 1:65){
  jpeg(paste0(names(p)[i], ".jpeg"), width=408, height=319)
  plot(p[[i]])
  dev.off()
}

setwd("..")

dir.create("cluster_plots")
dir.create("cluster_plots/sep_plots")
setwd("cluster_plots/sep_plots")

samps=c(1,5,6,9,10,20,25,50)

for(cl in unique(seuInt$cluster)){
  col_sp_clust=rep("gray80", length(unique(seuInt$cluster)))
  col_sp_clust[as.numeric(cl)]=cols_cluster[as.numeric(cl)]
  p_clust=SpaPlot(seuInt, batch=samps, point_size=1, cols=col_sp_clust, combine=F)
  
  p_clust=lapply(1:8, function(pl){p_clust[[pl]]+coord_fixed()+
      theme(legend.position="none")+ggtitle(metadata$sample[samps[pl]])})
  
  p=((p_clust[[1]] | p_clust[[2]] | p_clust[[3]] | p_clust[[4]]) /
       (p_clust[[5]] | p_clust[[6]] | p_clust[[7]] | p_clust[[8]]) & NoLegend()) +
    plot_annotation(title=paste0("Cluster ", cl)) &
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  jpeg(paste0("Cluster_", cl, ".jpeg"), width=600, height=382)
  plot(p)
  dev.off()
}

setwd("..")

pdf("clust_plots.pdf")

for(cl in sort(unique(seuInt$cluster))){
  col_sp_clust=rep("gray80", length(unique(seuInt$cluster)))
  col_sp_clust[as.numeric(cl)]=cols_cluster[as.numeric(cl)]
  p_clust=SpaPlot(seuInt, batch=samps, point_size=1, cols=col_sp_clust, combine=F)
  
  p_clust=lapply(1:8, function(pl){p_clust[[pl]]+coord_fixed()+
      theme(legend.position="none")+ggtitle(metadata$sample[samps[pl]])})
  
  p=((p_clust[[1]] | p_clust[[2]] | p_clust[[3]] | p_clust[[4]]) /
       (p_clust[[5]] | p_clust[[6]] | p_clust[[7]] | p_clust[[8]]) & NoLegend()) +
    plot_annotation(title=paste0("Cluster ", cl)) &
    theme(plot.title=element_text(hjust=0.5, face="bold"))
  
  plot(p)
}

dev.off()
