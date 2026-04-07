library(tidyverse)
library(glue)
library(pbmcapply)
setwd("/home/agelber/desp1/precast/prec_c25q25g3000/csrnp1_bam2/")
data_dir="/cndd2/agelber/hal/qc_aligned"
bam_fl="/home/agelber/desp1/precast/prec_c25q25g3000/csrnp1_bam2"

md=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/md_ser_all_new.rds") 

md=split(md, md$anot)

names(md)=make.names(names(md))

md=lapply(md, function(x) {
  split(x, x$sample)
})



dir.create("split_bams")
setwd("split_bams/")

for(clst in names(md)) {
  dir.create(clst)}

for (clst in names(md)) {
  
  
  
  for(smp in names(md[[clst]])) {
    
    data.table::fwrite(md[[clst]][[smp]] %>% dplyr::select(cell),
                       glue("{clst}/{smp}_barcodes.tsv"), 
                       sep = "\t", quote = F, 
                       col.names = F,
                       row.names = F)
    
    
    
  }
  
}

for (clst in names(md)) {
  
  
  invs=pbmclapply(names(md[[clst]]), function(smp){
    system(glue("/home/AD/agelber/desp1/subset-bam_linux ",
                "-b {bam_fl}/bams/{smp}.bam -c {clst}/{smp}_barcodes.tsv ",
                "-o {clst}/{smp}.bam --cores 20"))
    
    system(glue("rm {clst}/{smp}_barcodes.tsv"))
    
  })
  
  
}
