library(tidyverse)
library(Matrix)
library(patchwork)

setwd("~/desp1/precast/prec_c25q25g3000/figures/figure_1_plots/")

plt_objs=readRDS("~/desp1/shiny_app/plot_spatial/plt_objs.rds")

gns=c("Drd1",  "Trbc2", "Cux2","C1ql2","Agt", "Slc17a7","Neurod6", "Lamp5")

cols=pals::glasbey(6)[4]
to_pl=plt_objs$`042-B_WT_14 months_F`


plts=lapply(1:length(gns), function(i){
  coords=to_pl$coords
  coords[["gene"]]=as.numeric(to_pl$data[gns[i],])
  ggplot(coords)+geom_point(aes(x=col, y=-row, col=gene), size=0.01)+theme_void()+
    scale_color_gradient(name="", low="gray90", high=cols,
                         limits=c(0,
                                  quantile(coords$gene, .99)),
                         oob=scales::squish)+
    ggtitle(gns[i])+coord_fixed()+
    scale_y_reverse()+scale_x_reverse()+
    theme(plot.title=element_text(size = 16, hjust = 0.5),
          legend.title=element_text(size = 14),
          legend.text = element_text(size = 12))
  
})

p=wrap_plots(plts, ncol=3)
pdf("expr_plt.pdf",height = 4, width = 3)
plot(p)
dev.off()

seuInt=readRDS("~/desp1/precast/prec_c25q25g3000/seuInt_qc.rds")

md=seuInt@meta.data %>% dplyr::filter(sample=="042-B")

coords=to_pl$coords
img=to_pl$img
img=img[dim(img)[1]:1,dim(img)[2]:1,]
img=wrap_elements(grid::rasterGrob(as.raster(img)))

coords=coords %>% rownames_to_column("cell") %>%
  dplyr::filter(cell %in% md$cell)
coords=left_join(coords, md %>% select(cell, anot))

cols_cluster = c(ggsci::pal_material("green")(7)[2:7],
                 ggsci::pal_material("blue")(7)[2:7],
                 ggsci::pal_material("orange")(4)[2:4],
                 
                 ggsci::pal_material("red")(6)[2:6],
                 ggsci::pal_material("deep-purple")(4)[2:4]
)
names(cols_cluster)=readRDS("../cluster_order.rds")

hne_clusts=ggplot(coords)+
  geom_point(aes(x=-col, y=row, col=anot), size=2)+
  theme_void()+
  scale_color_manual(name="Precast Cluster",values = cols_cluster)+
  theme(legend.position = "none")+
  coord_fixed()+
  inset_element(img, left = 0, bottom = 0, right = 1, top = 1, on_top = F)

pdf("hne.pdf", paper = "a4r", compress = F)
plot(hne_clusts)
dev.off()
