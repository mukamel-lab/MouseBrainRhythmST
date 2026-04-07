library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(pathview)
library(enrichplot)

glist=readRDS("~/desp1/precast/prec_c25q25g3000/figures/fig3_app_wt_harm/dg_app_wt_glist.rds")
glist=split(glist, glist$class)
app_gs=mapIds(org.Mm.eg.db,glist$`APP23-TG`$gene,  'ENTREZID','SYMBOL')
app_gs=mapIds(org.Mm.egCHR,app_gs,  'chromosome','gene_id')

res_sig_app = readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/APP_harmonic/par_im/res_sig.rds")
res_sig_wt=readRDS("~/desp1/precast/prec_c25q25g3000/deseq_rhth/WT_harmonic/par_im/res_sig.rds")
res_sig_wt=split(res_sig_wt, res_sig_wt$cluster)
res_sig_app=split(res_sig_app, res_sig_app$cluster)

app_gs=lapply(names(res_sig_app), function(clt){
  mapIds(org.Mm.eg.db,
         setdiff(res_sig_app[[clt]]$gene, res_sig_wt[[clt]]$gene),
          'ENTREZID','SYMBOL')
  
})
names(app_gs)=names(res_sig_app)
kks=lapply(app_gs, function(gs) {
kk <- enrichKEGG(gene= gs,
organism     = 'mmu',
pvalueCutoff = 0.2)
})

res_kegg=lapply(kks, function(kk1){
  kk1@result
})
res_kegg=bind_rows(res_kegg, .id="cluster")
head(kk)
browseKEGG(kk, 'mmu04710')


mmu04068_dg <- pathview(gene.data  = app_gs$`Dentate Gyrus-sg`,
                     pathway.id = "mmu04068",
                     species    = "mmu"
                     )

mmu04068_dg <- pathview(gene.data  = gs,
                        pathway.id = "mmu04068",
                        species    = "mmu"
)

mmu04068_dg <- pathview(gene.data  = dg,
                        pathway.id = "mmu04010",
                        species    = "mmu"
)
