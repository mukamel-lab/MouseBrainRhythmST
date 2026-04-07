library(tidyverse)
library(glue)

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/csrnp1_bam")

md=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/md_ser_all_new.rds") %>%
  dplyr::filter(anot=="Cortex Layer 2/3")

md=split(md, md$sample)


dir.create("bams_ctx_l23")

setwd("bams")

for(smp in list.files()) {
  
  system(glue("samtools index {smp}"))
}

setwd("..")

dir.create("samp_barcodes")
setwd("samp_barcodes/")

for(nm in names(md)) {
  
  data.table::fwrite(md[[nm]] %>% dplyr::select(cell), paste0(nm, ".tsv"), sep = "\t", quote = F, col.names = F, row.names = F)
  
}

setwd("..")

for(nm in names(md)) {
  
system(glue("/home/AD/agelber/desp1/subset-bam_linux -b bams/{nm}.bam -c samp_barcodes/{nm}.tsv -o bams_ctx_l23/{nm}.bam"))

  
  }


meta.data=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/meta.data_filt.rds")

meta.data$grp=paste0(gsub(" months", "mo", meta.data$age), "_",
                     meta.data$genotype, "_",
                     meta.data$time)

meta.data=split(meta.data, meta.data$grp)

dir.create("merged_bams")
setwd("merged_bams/")
for (nm in names(meta.data)) {
  bams1=paste0("../bams_ctx_l23/",meta.data[[nm]]$sample,'.bam', collapse = " ")
  system(glue("samtools merge {nm}.bam {bams1}"))  
}

for (nm in names(meta.data)) {
  system(glue("samtools index {nm}.bam"))  
}
