library(glue)


bam_fl="/home/agelber/desp1/precast/prec_c25q25g3000/csrnp1_bam2"
data_dir="/cndd2/agelber/hal/qc_aligned"

for(fl in list.files(data_dir)){
  system(glue("samtools view  -b /cndd2/agelber/hal/qc_aligned/{fl}/outs/possorted_genome_bam.bam chr9:118000000-120000000 > {bam_fl}/bams/{fl}.bam"))
}

for(fl in list.files(data_dir)){
  system(glue("samtools index {bam_fl}/bams/{fl}.bam"))
}

