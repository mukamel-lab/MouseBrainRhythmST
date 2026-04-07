library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)                 
OrgDb = org.Mm.eg.db


res_interaction = readRDS("~/desp1/precast/precast_final_with_ros_caud/analysis2/rhythmicity_interaction/res_interaction.rds")
res_p=lapply(names(res_interaction), function(cl) {
  lapply(names(res_interaction[[cl]]), function(cmp) {
    y=res_interaction[[cl]][[cmp]]
    r1=substr(cmp,1,1)
    r2= substr(cmp,2,2)
    z=data.frame(gene=y$gene, 
                 cluster=cl,
                 rhythmicity_interaction_fdr=y$fdr,
                 region_1=r1, 
                 region_2=r2,
                 rhythmicity_fdr_region_1=y[[r1]],
                 rhythmicity_fdr_region_2=y[[r2]],
                 phi_region_1=y[[glue("phi_{r1}")]],
                 phi_region_2=y[[glue("phi_{r2}")]],
                 phi__hr_region_1=y[[glue("phi_hr_{r1}")]],
                 phi_hr_region_2=y[[glue("phi_hr_{r2}")]], 
                 amp_region_1=y[[glue("amp_{r1}")]],
                 amp_region_2=y[[glue("amp_{r2}")]]
                 
    )
  })%>% bind_rows()
})%>% bind_rows() 

res_p$region_1=plyr::mapvalues(res_p$region_1, from=c("R", "M"), to=c("rostral", "intermediate"))
res_p$region_2=plyr::mapvalues(res_p$region_2, from=c("M", "C"), to=c( "intermediate", "caudal"))
res_p$cntrst=paste0(res_p$region_1, "_vs_", res_p$region_2)


z=res_p %>% dplyr::filter(rhythmicity_fdr_region_1<0.1|rhythmicity_fdr_region_2<0.1)

z1=z %>% dplyr::filter(rhythmicity_interaction_fdr<0.1, region_1=="rostral",region_2=="caudal", 
                       amp_region_1>amp_region_2)


gene_symbols = unique(z1$gene)      
all_symbols  = unique(z$gene)     

fg_map  = bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
uni_map = bitr(all_symbols,  fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)

fg_entrez  = unique(fg_map$ENTREZID)
uni_entrez = unique(uni_map$ENTREZID)


ek = enrichKEGG(
  gene         = fg_entrez,
  organism     = "mmu",       # "mmu" for mouse, "hsa" for human
  universe     = uni_entrez,
  pvalueCutoff = 0.19,
  qvalueCutoff = .2
)


