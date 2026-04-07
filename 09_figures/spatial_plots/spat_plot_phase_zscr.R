library(tidyverse)

data_dir="/cndd2/agelber/hal/qc_aligned"

meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>%
mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>%
dplyr::filter(sample %in% list.files(data_dir)) %>% group_by(genotype) %>% tally

circ_var=function(x) {
 
    mean_sin = mean(sin(x))
    mean_cos = mean(cos(x))
    1 - sqrt(mean_sin^2 + mean_cos^2)
  
}

phi_dif= function(x,y) {
  atan2(sin(x-y),cos(x-y))}




res_sig <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/both_no_geno_term/par_im/res_sig.rds") %>%
  dplyr::filter(padj<0.01) %>% group_by(gene) %>% dplyr::filter(n()>=15)
shared_g=unique(res_sig$gene)
############################
#wts
##############################
wts=local({
coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/WT_harmonic/par_im/coefs.rds")

nms=names(coefs)

coefs=lapply(names(coefs), function(clust) {
  coefs[[clust]] %>% filter(gene %in% shared_g)   %>%
    mutate(phi=atan2(t_s,t_c)) %>% na.omit()
  
}) %>% setNames(nms) %>% bind_rows(.id="cluster")



coefs=coefs %>% group_by(gene) %>%
  mutate(zscr=phi_dif(phi, atan2(mean(t_s), mean(t_c)))/sqrt(circ_var(phi))) %>%
  na.omit() %>% ungroup() %>% group_by(cluster) %>% summarise(mean_zsq=mean(zscr^2)) 
})


############################
#all clusts sig APP
##############################

app=local({

coefs <- readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhyth2/APP_harmonic/par_im/coefs.rds")


nms=names(coefs)

coefs=lapply(names(coefs), function(clust) {
  coefs[[clust]] %>% filter(gene %in% shared_g)   %>%
    mutate(phi=atan2(t_s,t_c)) %>% na.omit()
  
}) %>% setNames(nms) %>% bind_rows(.id="cluster")



coefs=coefs %>% group_by(gene) %>%
  mutate(zscr=phi_dif(phi, atan2(mean(t_s), mean(t_c)))/sqrt(circ_var(phi))) %>%
  na.omit() %>% ungroup() %>% group_by(cluster) %>% summarise(mean_zsq=mean(zscr^2)) 


})


dta=readRDS("~/desp1/precast/prec_c25q25g3000/figures/spat_plots_generic/spat_data_030B.rds")
dta$cluster=as.character(dta$anot)

dta_wt=left_join(dta, wts)

dta_app=left_join(dta,app)
cluster_order= readRDS("~/desp1/precast/prec_c25q25g3000/figures/cluster_order.rds")
lmts=rbind(app, wts) %>%
  dplyr::filter(!(cluster %in% cluster_order[c(19:23,12)] )
)

lmts=max(lmts$mean_zsq)

dta_both=bind_rows(dta_wt %>% mutate(genotype="WT"), dta_app %>% mutate(genotype="APP23"))

plt=ggplot(dta_both %>%
             mutate(mean_zsq=ifelse(cluster %in% cluster_order[c(19:21,23,12)], NA, mean_zsq )))+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=mean_zsq), size=2)+
  theme_void(base_size = 16)+
  scale_color_scico(palette = 'oslo', name="Mean Squared\nCircular Z-Score",
                    oob=scales::squish) +
  ggtitle("Spatial patterns of phase variance")+coord_fixed()+
  theme(plot.title=element_text(size = 16, hjust = 0.5),
        legend.title=element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(.9, "lines")
  )+facet_wrap(~genotype)

plt

dta_both_dif=bind_rows(wts %>% mutate(genotype="WT"), app %>% mutate(genotype="APP23"))  %>% 
  ungroup() %>% group_by(cluster) %>%
  arrange(genotype) %>% 
  summarise(mean_zsq=mean_zsq[1]-mean_zsq[2])



dta_both_dif=left_join(dta, dta_both_dif)

lmts_dif= dta_both_dif%>%
  dplyr::filter(!(cluster %in% cluster_order[c(19:21,23,12)] )
  )
lmts_dif=max(abs(lmts_dif$mean_zsq))

plt_dif=ggplot(dta_both_dif %>%
             mutate(mean_zsq=ifelse(cluster %in% cluster_order[c(19:21,23,12)], NA, mean_zsq )))+
  geom_point(aes(y=-1*imagerow,x=imagecol, col=mean_zsq), size=2)+
  theme_void(base_size = 16)+
  scale_color_scico(palette = 'vik', name="Difference in\nMean Squared\nCircular Z-Score",midpoint = 0, 
                    limits=c(-1*lmts_dif,lmts_dif),
                    oob=scales::squish) +
  ggtitle("Spatial patterns of phase variance")+coord_fixed()+
  theme(plot.title=element_text(size = 16, hjust = 0.5),
        legend.title=element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(.9, "lines")
  )

plt_dif
