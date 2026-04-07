library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pathview)
library(enrichplot)
library(ReactomePA)

cluster_color <- readRDS("~/desp1/precast/prec_c25q25g3000/objects/cluster_color.rds")

res_sig= vroom::vroom("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_rhth_int/for_share/full/results.csv.gz") %>%
  dplyr::filter(padj<0.05)

ctx=res_sig %>% dplyr::filter(grepl("Cortex L", cluster))
ctx=res_sig %>% dplyr::filter(cluster %in% names(cluster_color)[1:6])

gs=mapIds(org.Mm.eg.db,unique(ctx$gene),  'ENTREZID','SYMBOL')

gs_hum=mapIds(org.Hs.eg.db,str_to_upper(unique(ctx$gene)),  'ENTREZID','SYMBOL')
kk= enrichKEGG(gene= gs,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)
edo2 <- pairwise_termsim(kk)
p1 <- emapplot(edo2 )

kk_hum= enrichKEGG(gene= gs_hum,
               organism = "hsa"  ,
               pvalueCutoff = 0.05)

kk_hum_all=enrichKEGG(gene= app_gs$`Cortex Layer 2/3`,
                      organism     = 'mmu',
                      pvalueCutoff = 0.5)

edo2 <- pairwise_termsim(kk_hum_all)
p1 <- emapplot(edo2 )
edo2 = enrichGO(cgs, OrgDb =org.Mm.eg.db ,
                qvalueCutoff = 0.05,ont="ALL")
go_hum=
edo = enrichDGN(gs)

edo = enrichGO(gs, OrgDb =org.Mm.eg.db ,qvalueCutoff = 0.1)
dotplot(edo2, showCategory=30) 
+ ggtitle("dotplot for ORA")

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

app_gs_hum=lapply(names(res_sig_app), function(clt){
  mapIds(org.Hs.eg.db,
         str_to_upper(setdiff(res_sig_app[[clt]]$gene, res_sig_wt[[clt]]$gene)),
         'ENTREZID','SYMBOL')
  
})
names(app_gs_hum)=names(res_sig_app)
ctx_appgs=unique(unlist(app_gs[names(cluster_color)[1:5]]))
ctx_appgs_hum=unique(unlist(app_gs_hum[names(cluster_color)[1:5]]))

edo2 = enrichGO(ctx_appgs, OrgDb =org.Mm.eg.db ,
                qvalueCutoff = 0.05,ont="ALL")
dotplot(edo2, showCategory=15) + ggtitle("GO for CTX Gain of Rhythmicity")

edo2 = enrichDO(ctx_appgs_hum ,
                qvalueCutoff = 0.05)
dotplot(edo2, showCategory=15) + ggtitle("GO for CTX Gain of Rhythmicity")

edo2 = enrichGO(ctx_appgs, OrgDb =org.Mm.eg.db ,
                qvalueCutoff = 0.05,ont="ALL")
dotplot(edo2, showCategory=15) + ggtitle("GO for CTX Gain of Rhythmicity")

edo2 <- pairwise_termsim(edo2)
p1 <- emapplot(edo2, )

upsetplot(kk)

edo = enrichDGN(app_gs_hum$`Cortex Layer 2/3`)

dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")

edo2 = enrichGO(app_gs$`Cortex Layer 2/3`, OrgDb =org.Mm.eg.db ,
                qvalueCutoff = 0.1,ont="ALL")
dotplot(edo2, showCategory=20) + ggtitle("GO for CTX L2/3")







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
browseKEGG(kk, 'mmu04723')


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

mmu04068_dg <- pathview(gene.data  = gs,
                        pathway.id = "mmu04068",
                        species    = "mmu"
)
