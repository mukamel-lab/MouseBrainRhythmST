library(tidyverse)
library(glue)

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/csrnp1_bam2")

meta.data=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/meta.data_filt.rds")

meta.data$grp=paste0(gsub(" months", "mo", meta.data$age), "_",
                     meta.data$genotype, "_",
                     meta.data$time)

meta.data=split(meta.data, meta.data$grp)


setwd("merged_bams/")
clsts=list.files("../split_bams")
for(clst in clsts) {
  dir.create(clst)
for (nm in names(meta.data)) {
  bams1=paste0("../split_bams/",clst,"/",meta.data[[nm]]$sample,'.bam', collapse = " ")
  
  system(glue("samtools merge {clst}/{nm}.bam {bams1}"))  
  system(glue("samtools index {clst}/{nm}.bam"))  
}


}
