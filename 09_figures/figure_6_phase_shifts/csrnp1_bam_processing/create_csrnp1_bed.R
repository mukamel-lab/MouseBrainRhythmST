library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)

data_dir="/cndd2/agelber/hal/qc_aligned"

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/csrnp1_bam")

txdb = makeTxDbFromGFF("/home/AD/agelber/space/genes.gtf", format = "gtf")
csrnp1=as.data.frame(import.gff("/home/AD/agelber/space/genes.gtf", format = "gtf"))
csrnp1= csrnp1 %>% filter(grepl("Csrnp1", gene_name))

genes = genes(txdb, filter = list(GENEID  = unique(csrnp1$gene_id)))

export.bed(genes, "csrnp1.bed")

gene_ranges = import.bed("csrnp1.bed")

bam_file ="/cndd2/agelber/hal/qc_aligned/030-B/outs/possorted_genome_bam.bam"
#bam = open(BamFile(bam_file))
reads_for_gene = scanBam(bam_file,index = bam_file, which = gene_ranges)

filtered_bam =  "test.bam"
filtered_bam_file = BamFile(filtered_bam, open = "w", index = TRUE)
writeBamFile(reads_for_gene, filtered_bam_file)
close(filtered_bam_file)
close(bam)
