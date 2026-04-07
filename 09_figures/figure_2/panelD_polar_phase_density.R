library(tidyverse)
library(patchwork)

cluster_color <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

compute_dens_polar=function(pho) {

    phi.df=tibble(phi=rep((1:1000-1)/500*pi,1))
    phi.df$hour=phi.df$phi/pi*12
    phi.df$count=phi.df$phi
    kapp=20
    
    for(s in 1:length(phi.df$phi)){
      phi.df$count[s]=sum(exp(kapp*cos(phi.df$phi[s]-as.numeric(pho))))
    }
    phi.df$dens=phi.df$count/sum(phi.df$count)/phi.df$phi[2] #normalize in the linear sense
    phi.df$dens=phi.df$count/sqrt(sum(phi.df$count^2)*phi.df$phi[2]/2)#to normalize the radar plot integral to one

   return(phi.df)
  
  }


############################
 # WT:  all clusts sig 
##############################

res_sig <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res_sig.rds")
res_sig=split(res_sig, res_sig$cluster)
coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/coefs.rds")


phi_dens_wt=lapply(names(coefs), function(clust) {
  coefs_clust=coefs[[clust]] %>% filter(gene %in% res_sig[[clust]]$gene)
  compute_dens_polar(coefs_clust$phi)
})
names(phi_dens_wt)=names(coefs)
phi_dens_wt=bind_rows(phi_dens_wt, .id="cluster")

plt=ggplot(phi_dens_wt %>% mutate(cluster=factor(cluster, levels=names(cluster_color))))+geom_line(aes(x=hour,y=dens, color=cluster))+
  scale_color_manual(values= cluster_color, name="Cluster")+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
  theme_minimal()+labs(x="ZT", y="Density", title=paste("Phase Density CTX and DG Sig Genes"))+
  theme(legend.position = NULL, text = element_text(size=12))  


phi_dens_wt$cluster[phi_dens_wt$cluster=="Olfactory Tubercle"]="Amygdala"
names(cluster_color)[names(cluster_color)=="Olfactory Tubercle"]="Amygdala"
phi_dens_wt$region=plyr::mapvalues(phi_dens_wt$cluster,names(cluster_color),
                                c(rep("Cortex", 6), rep("Hippocampus",6),
                                  rep("Olfactory Areas", 3), rep("Cerebral Nuclei", 5),
                                  rep("Non-neuronal", 3)))
phi_dens_wt$region=factor(phi_dens_wt$region, levels=c("Cortex","Hippocampus","Olfactory Areas",
                                                 "Cerebral Nuclei","Non-neuronal"))

phi_dens_wt$geno="NTG"
plt_wt=ggplot(phi_dens_wt %>% mutate(cluster=factor(cluster, levels=names(cluster_color))))+geom_line(aes(x=hour,y=dens, color=cluster))+
  scale_color_manual(values= cluster_color, name="Cluster")+coord_polar()+scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+ylim(0,NA)+
  theme_minimal()+labs(x="ZT", y="Density", title=paste("Phase Density Across Regions for Significantly Rhythmic Genes"))+
  theme(legend.position = NULL, text = element_text(size=16))+facet_wrap(~region)

plt_wt


############################
# APP:  all clusts sig 
##############################


res_sig <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/res_sig.rds")
res_sig=split(res_sig, res_sig$cluster)
coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/coefs.rds")


phi_dens_app=lapply(names(coefs), function(clust) {
  coefs_clust=coefs[[clust]] %>% filter(gene %in% res_sig[[clust]]$gene)
  compute_dens_polar(coefs_clust$phi)
})
names(phi_dens_app)=names(coefs)
phi_dens_app=bind_rows(phi_dens_app, .id="cluster")



phi_dens_app$cluster[phi_dens_app$cluster=="Olfactory Tubercle"]="Amygdala"
names(cluster_color)[names(cluster_color)=="Olfactory Tubercle"]="Amygdala"
phi_dens_app$region=plyr::mapvalues(phi_dens_app$cluster,names(cluster_color),
                                   c(rep("Cortex", 6), rep("Hippocampus",6),
                                     rep("Olfactory Areas", 3), rep("Cerebral Nuclei", 5),
                                     rep("Non-neuronal", 3)))
phi_dens_app$region=factor(phi_dens_app$region, levels=c("Cortex","Hippocampus","Olfactory Areas",
                                                       "Cerebral Nuclei","Non-neuronal"))

phi_dens_app$geno="APP23"
plt_app=ggplot(phi_dens_app %>% mutate(cluster=factor(cluster, levels=names(cluster_color))))+
  geom_line(aes(x=hour,y=dens, color=cluster))+
  scale_color_manual(values= cluster_color, name="Cluster")+
  coord_polar()+
  scale_x_continuous(breaks=seq(0, 24, by=4),expand=c(0,0), lim=c(0, 24))+
  ylim(0,NA)+
  theme_minimal()+
  labs(x="ZT", y="Density", title=paste("Phase Density Across Regions for Significantly Rhythmic Genes"))+
  theme(legend.position = NULL, text = element_text(size=16))+facet_wrap(~region)

plt_app

############################
# APP:  not sig in WT
##############################

res_sig_wt <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res_sig.rds")
res_sig_app <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/res_sig.rds")
coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/coefs.rds")

# Split by cluster for easier access
res_sig_wt = split(res_sig_wt, res_sig_wt$cluster)
res_sig_app = split(res_sig_app, res_sig_app$cluster)

# Filter APP genes that are not significant in WT
phi_dens_app_ns_wt <- lapply(names(res_sig_app), function(clust) {
  app_genes = res_sig_app[[clust]]$gene
  wt_genes = res_sig_wt[[clust]]$gene
  app_only_genes = setdiff(app_genes, wt_genes)  # Genes in APP but not in WT
  coefs_clust = coefs[[clust]] %>% filter(gene %in% app_only_genes)
  compute_dens_polar(coefs_clust$phi)
})
names(phi_dens_app_ns_wt) = names(coefs)
phi_dens_app_ns_wt = bind_rows(phi_dens_app_ns_wt, .id="cluster")

phi_dens_app_ns_wt$cluster[phi_dens_app_ns_wt$cluster=="Olfactory Tubercle"]="Amygdala"
names(cluster_color)[names(cluster_color)=="Olfactory Tubercle"]="Amygdala"
phi_dens_app_ns_wt$region=plyr::mapvalues(phi_dens_app_ns_wt$cluster,names(cluster_color),
                                    c(rep("Cortex", 6), rep("Hippocampus",6),
                                      rep("Olfactory Areas", 3), rep("Cerebral Nuclei", 5),
                                      rep("Non-neuronal", 3)))
phi_dens_app_ns_wt$region=factor(phi_dens_app_ns_wt$region, levels=c("Cortex","Hippocampus","Olfactory Areas",
                                                         "Cerebral Nuclei","Non-neuronal"))

phi_dens_app_ns_wt$geno="APP23"

# Plot
plt_app_ns_wt = ggplot(phi_dens_app_ns_wt %>% mutate(cluster=factor(cluster, levels=names(cluster_color)))) +
  geom_line(aes(x=hour, y=dens, color=cluster)) +
  scale_color_manual(values= cluster_color, name="Cluster") +
  coord_polar() +
  scale_x_continuous(breaks=seq(0, 24, by=4), expand=c(0, 0), lim=c(0, 24)) +
  ylim(0, NA) +
  theme_minimal() +
  labs(x="ZT", y="Density", title="Phase Density: APP Genes Not Significant in WT") +
  theme(legend.position = NULL, text = element_text(size=12))+facet_wrap(~region)


plt_app_ns_wt

############################
# WT:  not sig in APP
##############################

coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/coefs.rds")

phi_dens_wt_ns_app <- lapply(names(res_sig_wt), function(clust) {
  wt_genes = res_sig_wt[[clust]]$gene
  app_genes = res_sig_app[[clust]]$gene
  wt_only_genes = setdiff(wt_genes, app_genes)  # Genes in WT but not in APP
  coefs_clust = coefs[[clust]] %>% filter(gene %in% wt_only_genes)
  compute_dens_polar(coefs_clust$phi)
})
names(phi_dens_wt_ns_app) = names(coefs)
phi_dens_wt_ns_app = bind_rows(phi_dens_wt_ns_app, .id="cluster")

phi_dens_wt_ns_app$cluster[phi_dens_wt_ns_app$cluster=="Olfactory Tubercle"]="Amygdala"
names(cluster_color)[names(cluster_color)=="Olfactory Tubercle"]="Amygdala"
phi_dens_wt_ns_app$region=plyr::mapvalues(phi_dens_wt_ns_app$cluster,names(cluster_color),
                                          c(rep("Cortex", 6), rep("Hippocampus",6),
                                            rep("Olfactory Areas", 3), rep("Cerebral Nuclei", 5),
                                            rep("Non-neuronal", 3)))
phi_dens_wt_ns_app$region=factor(phi_dens_wt_ns_app$region, levels=c("Cortex","Hippocampus","Olfactory Areas",
                                                                     "Cerebral Nuclei","Non-neuronal"))

phi_dens_wt_ns_app$geno="NTG"


# Plot
plt_wt_ns_app = ggplot(phi_dens_wt_ns_app %>% mutate(cluster=factor(cluster, levels=names(cluster_color)))) +
  geom_line(aes(x=hour, y=dens, color=cluster)) +
  scale_color_manual(values= cluster_color, name="Cluster") +
  coord_polar() +
  scale_x_continuous(breaks=seq(0, 24, by=4), expand=c(0, 0), lim=c(0, 24)) +
  ylim(0, NA) +
  theme_minimal() +
  labs(x="ZT", y="Density", title="Phase Density: WT Genes Not Significant in APP") +
  theme(legend.position = NULL, text = element_text(size=12))+facet_wrap(~region)

plt_wt_ns_app

