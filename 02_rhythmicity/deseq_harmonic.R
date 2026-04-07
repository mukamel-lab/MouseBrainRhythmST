library(pbmcapply)
library(DESeq2)
library(tidyverse)
library(glmGamPoi)
library(customfuncs)


setwd("~/desp1/precast/prec_c25q25g3000/")
data_dir="/cndd2/agelber/hal/qc_aligned"

  agg_exp = readRDS("objects/agg_c.rds")
  
meta.data = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
    mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
    dplyr::filter(sample %in% list.files(data_dir)) 
  
meta.data$t_c=as.numeric(gsub("ZT", "", meta.data$time))
meta.data$t_s=as.numeric(gsub("ZT", "", meta.data$time))
meta.data= meta.data%>% mutate(t_s=round(sin(2*pi*(t_s)/24)),
                               t_c=round(cos(2*pi*(t_c)/24)))


deseqs=pbmclapply(agg_exp, function(cl){
  # full model specified in design
  # reduced model specified in reduced
  try({
    suppressWarnings({
      counts1=cl
    
      md=meta.data %>% dplyr::rename(id=sample)
      maxs=apply(counts1 %>% dplyr::select(-gene) , 1, max)
      counts2=counts1[maxs>10, ]
      cs=colSums(counts2 %>% select(-gene))
      cs=names(cs)[cs>2e5]
      ngens=apply(counts2 %>% dplyr::select(-gene), 2, function(x){sum(x>0)})
      ngens=ngens/nrow(counts2)
      ngens=names(ngens)[ngens>0.9]
      
      md=md %>% filt(id %in% intersect(cs, ngens), genotype=="WT")
      
      DESeqDataSetFromMatrix(countData=counts2[,c("gene", md$id)] %>% 
                               filt(!(gene %in% c("humanAPP", "Thy1") )), 
                             colData=md, 
                             design=~1+t_s+t_c,
                             tidy = T) %>%
        estimateSizeFactors(type="poscounts") %>%
        estimateDispersions(fitType="glmGamPoi") %>%
        DESeq(test="LRT", reduced=~1)
      
      
    })
  })
}, mc.allow.recursive = T, ignore.interactive = T)


deseqs=deseqs[which(sapply(deseqs, function(x) {is(x, "DESeqDataSet")}))]

res=lapply(deseqs, function(x) {results(x, tidy = T)}) %>% 
  bind_rows(.id="cluster") %>%
  dplyr::rename(gene=row)
res_sig= res %>% subset(padj<0.05)

ncounts=lapply(deseqs, 
               function(x) { counts(x, normalized=T) %>%
                   as.data.frame() %>% 
                   rownames_to_column("gene") %>%
                   pivot_longer(-gene, names_to = "sample", values_to = "l2expr") %>%
                   left_join(meta.data %>%
                               dplyr::select(-t_c, -t_s))})


ncounts=lapply(deseqs, 
               function(x) {counts(x, normalized=T) %>%
                   as.data.frame() %>% 
                   rownames_to_column("gene") %>%
                   pivot_longer(-gene, names_to = "sample", values_to = "norm_counts") %>%
                   left_join(meta.data %>%
                               dplyr::select(-t_c, -t_s)) %>%
                   mutate(l2expr=log2(norm_counts+1))
                 
                 })

coefs=lapply(deseqs, 
               function(x) { coefficients(x) %>%
                   as.data.frame() %>% 
                   rownames_to_column("gene") %>%
                   dplyr::rename(a=t_c, b=t_s) %>%
                   mutate(phi=atan2(b,a) %% (2*pi), phi_hr=phi*12/pi, amp=sqrt(a^2+b^2))})


coefs_counts=lapply(names(deseqs), function(cl){
  x=ncounts[[cl]]
  y=coefs[[cl]]
  z=left_join(x,y, by="gene")
  z
})

names(coefs_counts)=names(deseqs)

setwd("deseq_rhth")
par_saveim("WT_harmonic_basic")

setwd("WT_harmonic_basic")


coefs_counts=bind_rows(coefs_counts,.id="cluster")

all_data=left_join(coefs_counts, res, by=c("gene", "cluster"))

all_data$rel_amp=all_data$amp/log2(all_data$baseMean+1)

data.table::fwrite(all_data %>% na.omit(), "all_results.csv.gz", row.names = F)

dir.create("par_im")

system("mv *rds par_im")