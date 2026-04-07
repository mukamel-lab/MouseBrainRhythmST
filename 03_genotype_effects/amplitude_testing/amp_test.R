library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)
library(patchwork)
library(MASS)

setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

agg_exp = readRDS("objects/agg_c.rds")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")


meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir)) 

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))
# deseqs = pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl_n) {
#   # full model specified in design
#   # reduced model specified in reduced
#   try({
#     suppressWarnings({
#       
#       counts1 = agg_exp[[cl_n]]
#       
#       md = meta.data %>% dplyr::rename(id = sample) %>%
#         dplyr::filter(age == "7 months")
#       
#       counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
#       
#       maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
#       counts2 = counts1[maxs > 4, ]
#       
#       exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
#       exp_gns = exp_gns[exp_gns > 0.8]
#       
#       md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
#       
#      DESeqDataSetFromMatrix(
#         countData = counts2[, c("gene", md$id)]  ,
#         colData = md,
#         design =  ~ genotype + sex + t_s + t_c,
#         tidy = T
#       ) %>%
#         estimateSizeFactors(type = "ratio") %>%
#         estimateDispersions(fitType = "local") %>%
#         DESeq(test = "LRT",
#               reduced =  ~ genotype+sex,
#               fitType = "local")})})})
# 
# saveRDS(deseqs,"~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/amplitude_testing/deseqs_young.rds" )
deseqs=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/amplitude_testing/deseqs_young.rds" )
setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new")

cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)

res_wt_aw = readRDS("young/res_y_aw.rds") %>%
  dplyr::filter( padj<0.1) %>%  group_by(cluster, gene) %>% 
  dplyr::filter(n()==2) 

res_wt_aw=split(res_wt_aw, res_wt_aw$cluster) 

coefs = readRDS("young/coefs_y.rds") 

coefs=split(coefs, coefs$cluster)

coefs=lapply(names(res_wt_aw), function(clst) {
  cf=coefs[[clst]] %>% 
    dplyr::filter(gene %in% res_wt_aw[[clst]]$gene )
  cf
}) %>% bind_rows()


coefs=coefs %>%
  mutate(rel_amp_dif=app_amp-wt_amp, 
         phi_dif=atan2(sin(app_phi-wt_phi),cos(app_phi-wt_phi)) ,
phi_dif_hr=(12*phi_dif/pi+24)%%24)

coefs_sim=coefs %>% dplyr::filter(abs(phi_dif_hr)<=2)
coefs_sim=split(coefs_sim, coefs_sim$cluster)

genes_by_clust=lapply(coefs_sim, function(x) {x$gene})


coefs_joint=lapply(deseqs, function(x) {
   coefficients(x) %>% as.data.frame() %>% rownames_to_column("gene")%>% 
    mutate(phi =  atan2(t_s, t_c) %% (2 * pi)) 
})


tasks <- data.frame(
  cluster = rep(names(genes_by_clust), times = lengths(genes_by_clust)),
  gene    = unlist(genes_by_clust, use.names = FALSE)
)
tasks=split(tasks,tasks$cluster)

counts_all=lapply(deseqs, counts)
res_list=vector("list", length(tasks))
names(res_list)=names(tasks)
ij=1

for(tsk in tasks){

    try({
      cluster=tsk$cluster[1]
      cnts=counts_all[[cluster]]
      cdta=as.data.frame(deseqs[[cluster]]@colData)
      cfj=coefs_joint[[cluster]] 
      phi_vec <- setNames(cfj$phi,
                        cfj$gene)
      thetas=dispersions(deseqs[[cluster]])
      thetas=1/thetas
      names(thetas)=rownames(cnts)
     
      
      
 res_c=lapply(1:nrow(tsk), function(i){

      gene    <- tsk$gene[i]
    # get phase for this gene
   
    phi     <- phi_vec[gene]
    
    # build per‐gene data.frame
    df <- tibble(
      id     = names(cnts[gene,]),
      counts = cnts[gene,]
    ) %>%
      left_join(cdta %>% as.data.frame(), by = "id") %>%
      mutate(
        time     = pi * as.numeric(sub("ZT", "", time)) / 12,
        harmonic = cos(time - phi)
      )
    
    fit_full <- glm.nb(
      counts ~ genotype * harmonic + offset(log(sizeFactor)),
      data    = df,
      init.theta = thetas[gene]
    )
    


    fit_red  <-  update(fit_full, . ~ . - genotype:harmonic)

    
    # extract LRT & p‐value
    LL1  <- as.numeric(logLik(fit_full))
    LL2  <- as.numeric(logLik(fit_red))
    LRT  <- 2 * (LL1 - LL2)
    pval <- pchisq(LRT, df = 1, lower.tail = FALSE)
    
    # return one‐row tibble
    data.frame(
      cluster    = cluster,
      gene       = gene,
      full_LL    = LL1,
      reduced_LL = LL2,
      LRT        = LRT,
      p.value    = pval
    )
    
    }) %>% bind_rows() %>% ungroup() %>%
   mutate(padj=p.adjust(p.value, method="BH"))
     res_c=left_join(res_c, 
                     coefs_sim[[cluster]] %>%
                       dplyr::select(gene,rel_amp_dif) %>% dplyr::rename(amp_dif=rel_amp_dif))
     res_c=left_join(res_c, res_wt_aw[[cluster]] %>% dplyr::select(gene, padj, genotype) %>% 
                       pivot_wider(id_cols = "gene", names_from = "genotype", values_from = "padj", 
                                                                                               names_prefix ="rhth_" ))
       res_list[[cluster]]=res_c
    })
  print(ij)
  ij=ij+1
}

results_df <- bind_rows(res_list)
saveRDS(results_df, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/amplitude_testing/res_list_young2.rds")

View(results_df)
