library(tidyverse)
library(glue)
library(ComplexHeatmap)
library(magrittr)
library(scico)
library(viridis)

setwd("/home/agelber/desp1/precast/prec_c25q25g3000/")

res = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/deseq_reg_int2/regional_DRG_NTG.rds") %>% 
  bind_rows(.id="comparison")

fdr_cut=0.05

short_nms=readRDS("objects/short_names2.rds")

names(short_nms)[1]="L23"

jci=matrix(0, nrow=23, ncol=23)

rownames(jci)=names(short_nms)
colnames(jci)=names(short_nms)


for(i in 1:23) {
  for(j in 1:23) {
    
    if(i==j){
      jci[i,j]=NA
    } else if (i>j) {
      
      df=res %>% dplyr::filter((region_1== rownames(jci)[i] | region_2==rownames(jci)[i] ) &
                               (region_1== rownames(jci)[j] | region_2==rownames(jci)[j] )
                               )
      
      jci[i,j]=sum(df$rhythmic_in_r1*df$rhythmic_in_r2)/nrow(df)
    } else {
      
      df=res %>% dplyr::filter((region_1== rownames(jci)[i] | region_2==rownames(jci)[i] ) &
                                 (region_1== rownames(jci)[j] | region_2==rownames(jci)[j] )) 
      # df=res %>%
      #   dplyr::filter(rhythmic_in_r1, rhythmic_in_r2)
      jci[i,j]=log10(sum(df$fdr<fdr_cut, na.rm = T)+1)
      
    }
    
  }
}


##############################
logN_vals = jci[upper.tri(jci)]
jac_vals  = jci[lower.tri(jci)]

logN_max = max(logN_vals, na.rm = TRUE)
jac_max  = max(jac_vals,  na.rm = TRUE)
sev_pct=quantile(logN_vals,0.7)


col_fun1 = circlize::colorRamp2(
  seq(0, logN_max, length.out = 1000),
  viridis::viridis(1000)
)

col_fun2 = circlize::colorRamp2(
  seq(0, jac_max, length.out = 1000),
  viridis::magma(1000)
)

hm_wt = Heatmap(name="hm_wt",
                 jci, 
                 cluster_rows = FALSE, 
                 cluster_columns = FALSE, 
                 show_heatmap_legend = FALSE,
                 na_col = "white",
                 col = col_fun1, 
                 column_names_side = "bottom",
                 row_names_side = "left",
                column_names_gp = gpar(fontsize = 6, col = "black",fontfamily="ArialMT"),
                row_names_gp = gpar(fontsize = 6, col = "black",fontfamily="ArialMT"),
                
                 rect_gp = gpar(type = "none"),
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if (i > j) {
                     fill = col_fun2(jci[i, j])
                     grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                     if(jci[i,j]>=.2) {
                       grid.text( signif(jci[i,j],2), x, y, gp = gpar(fontsize = 6, col = "black",fontfamily="ArialMT" ))
                     }
                   } else {
                     fill = col_fun1(jci[i, j])
                     grid.rect(x, y, w, h, gp = gpar(fill = fill, col = NA))
                     if(!is.na(jci[i,j]) & jci[i,j]>sev_pct) {
                       grid.text( round(10^jci[i,j]), x, y, gp = gpar(fontsize = 6, col = "black",fontfamily="ArialMT" ))
                     }
                   }
                 }
)
# Draw the heatmap

legend1 = Legend(
  at = c(0, logN_max),
  col_fun = col_fun1,
  title = "Number of DRG (FDR<0.05)",
  labels_gp = gpar(fontsize = 7, col = "black",
                   fontfamily="ArialMT"),
  title_gp = gpar(fontsize = 7, col = "black",
                  fontfamily="ArialMT", fontface="bold"),
  title_position = "leftcenter-rot",
  labels = c("0", as.character(round(10^logN_max - 1)))
)

legend2 = Legend(
  at = c(0, jac_max),
  col_fun = col_fun2,
  labels_gp = gpar(fontsize = 7, col = "black",
                   fontfamily="ArialMT"),
  title_gp = gpar(fontsize = 7, col = "black",
                  fontfamily="ArialMT", fontface="bold"),
  title_position = "leftcenter-rot",
  title = "Jaccard index",
  labels = c("0", signif(jac_max, 2))
)


draw(hm_wt, padding = unit(c(0, 0, 0, 3.5), "cm"))

for(i in c(6,12,15,20)){
  
  decorate_heatmap_body("hm_wt", 
                        {grid.lines(c(i/23,i/23),
                                    c(0,
                                      1),
                                    gp = gpar(lty = 2, lwd = 1.5, 
                                              col="gray85"))
                        })
  
  decorate_heatmap_body("hm_wt", 
                        {grid.lines(c(0, 1),
                                    c((23-i)/23, 
                                      (23-i)/23),
                                    gp = gpar(lty = 2, lwd = 1.5, 
                                              col="gray85"))
                        })
}
  
  

pushViewport(viewport(x = 0.82, y = 0.65, width = 0.2,
                      height = 0.5, 
                      just = c("left", "center")))
draw(legend1)

popViewport()

pushViewport(viewport(x = 0.82, y = 0.35, width = 0.2,
                      height = 0.5,
                      just = c("left", "center")))
draw(legend2)
popViewport()







