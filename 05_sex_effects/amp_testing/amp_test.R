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
  dplyr::filter(sample %in% list.files(data_dir)) %>% dplyr::filter(genotype=="WT")

meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))
deseqs = pbmclapply(names(agg_exp) %>% setNames(.,.), function(cl_n) {
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      
      counts1 = agg_exp[[cl_n]]
      
      md = meta.data %>% dplyr::rename(id = sample)
      
      counts1 = counts1[, c("gene", intersect(colnames(counts1), md$id))]
      
      maxs = rowSums((counts1 %>% dplyr::select(-gene)) > 10)
      counts2 = counts1[maxs > 4, ]
      
      exp_gns = colSums(counts2 %>% dplyr::select(-gene) > 0) / nrow(counts2)
      exp_gns = exp_gns[exp_gns > 0.8]
      
      md = md %>% dplyr::filter(id %in% colnames(counts2), id %in% names(exp_gns))
      
     DESeqDataSetFromMatrix(
        countData = counts2[, c("gene", md$id)]  ,
        colData = md,
        design =  ~ age+sex + t_s + t_c,
        tidy = T
      ) %>%
        estimateSizeFactors(type = "ratio") %>%
        estimateDispersions(fitType = "local") %>%
        DESeq(test = "LRT",
              reduced =  ~ age + sex,
              fitType = "local")})})})



cluster_color = readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")
short_names = readRDS("~/desp1/precast/prec_c25q25g3000/objects/short_names2.rds")
names(cluster_color)=plyr::mapvalues(names(cluster_color), short_names, names(short_names))

regions=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_order2.rds")

regions=setNames(names(short_names), names(regions))
region_color=readRDS("~/desp1/precast/prec_c25q25g3000/objects/region_color.rds")
names(region_color)=names(regions)


coefs=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/interaction_results_readjusted.rds")
coefs=coefs%>% mutate(phi_dif= atan2(sin(f_phi-m_phi),cos(f_phi-m_phi)) ,phi_dif_hr=12*phi_dif/pi)

coefs_sim=coefs %>% dplyr::filter(abs(phi_dif_hr)<=1)
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
      counts ~ sex * harmonic + offset(log(sizeFactor)),
      data    = df,
      init.theta = thetas[gene]
    )
    
   

    fit_red  <-  update(fit_full, . ~ . - sex:harmonic)

    
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
       res_list[[cluster]]=res_c
    })

}

results_df <- bind_rows(res_list)
saveRDS(results_df, "~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/sex_int_geno_spec/amp_testing/res_df.rds")

View(results_df)
