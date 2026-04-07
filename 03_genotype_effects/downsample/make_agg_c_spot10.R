library(Matrix)
library(tidyverse)
library(customfuncs)
library(glue)
library(pbmcapply)
library(Matrix.utils)

setwd("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/geno_int_new/downsample")


metadata = readRDS("~/desp1/precast/prec_c25q25g3000/objects/md_ser_all_new.rds") %>% 
  group_by(anot,sample, time, genotype, age) %>% tally() %>% dplyr::filter(n>=35) %>%
  group_by(anot, time, genotype, age) %>% tally() %>% dplyr::filter(n>=2) %>% 
  group_by(anot, genotype, age) %>% tally() %>% dplyr::filter(n==4) %>% 
  group_by(anot) %>% tally()  %>% dplyr::filter(n==4) %>% 
  pull(anot) %>% unique()

raw_counts_list=readRDS("raw_counts_list.rds")

raw_counts_list=raw_counts_list[metadata]
for(i in 1:30) {
  set.seed(i)
agg_c=pbmclapply(raw_counts_list, function(smp){
  set.seed(i)
  md=smp$md %>% rownames_to_column("rn")
  md=md %>% group_by(sample) %>% dplyr::filter(n()>=30)%>%
    sample_n(25, replace =F)
  
  
  x=aggregate.Matrix(t(smp$counts[,md$rn]), groupings = md$sample, fun="sum") %>%
    as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("gene")
  
  return(x)})

saveRDS(agg_c,
        glue("~/desp1/precast/prec_c25q25g3000/deseq_rhth_int2/",
        "geno_int_new/downsample/spot25/agg_c_ds{i}.rds"))
}