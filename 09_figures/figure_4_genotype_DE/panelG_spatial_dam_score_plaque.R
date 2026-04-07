library(tidyverse)
library(scattermore)
library(scales)
logcpm_dam=readRDS("~/desp1/precast/prec_c25q25g3000/figures/fig4_app_basic/dam_scr/logcpm_dam.rds")
logcpm_dam=apply(logcpm_dam, 2, function(x) {log10((1e6*x/sum(x)+1))})
logcpm_dam=t(logcpm_dam)
logcpm_dam=logcpm_dam[which(!is.na(logcpm_dam[,2])),]

pc = prcomp(logcpm_dam,
             center = TRUE,
             scale. = TRUE)
md_ser_all = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/objects/md_ser_all.rds")

md_ser_all=md_ser_all[rownames(pc$x),]
md_ser_all=cbind(md_ser_all, pc$x)

coords=readRDS("~/desp1/precast/prec_c25q25g3000/figures/fig4_app_basic/dam_scr/coords_df.rds")
ggplot(md_ser_all)+geom_scattermore(aes(x=PC1,y=PC2, col=genotype))#+facet_wrap(~anot)
md_ser_all=left_join(md_ser_all,coords, by=c("sample","cell"))
#md_ser_all = readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/fig4_app_basic/dam_scr/md_pc_mat_cpm_only_dam.rds")
md_ser_all=md_ser_all %>% ungroup() %>% 
  mutate(dam_scr=-1*(PC1+PC2))
scale_lims=c(min(md_ser_all$dam_scr),max(md_ser_all$dam_scr))

scale_lims=c(quantile(md_ser_all$dam_scr,.001),quantile(md_ser_all$dam_scr,.999))

md_split=split(md_ser_all, md_ser_all$sample)

data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir))  %>% arrange(sex)

pdf("young_wt.pdf")
for(df in md_split[metadata$sample[metadata$age=="7 months" & metadata$genotype=="WT"]]) {
p=ggplot(df)+ geom_point(aes(x=img_x,y=img_y,col=dam_scr))+coord_fixed()+theme_void()+
  ggtitle(paste0(c(df$sample[1], df$genotype[1], df$age[1], df$sex[1]), collapse = "_"))+
  theme(plot.title = element_text(hjust=0.5))+
  scale_color_gradient2(name="DAM Score", 
                        low = "gray30", 
                        high ="red",
                        limits=c(scale_lims), 
                        oob=scales::squish)
plot(p)

}

dev.off()

pdf("young_app.pdf")
for(df in md_split[metadata$sample[metadata$age=="7 months" & metadata$genotype!="WT"]]) {
  p=ggplot(df)+ geom_point(aes(x=img_x,y=img_y,col=dam_scr))+coord_fixed()+theme_void()+
    ggtitle(paste0(c(df$sample[1], df$genotype[1], df$age[1], df$sex[1]), collapse = "_"))+
    theme(plot.title = element_text(hjust=0.5))+
    scale_color_gradient2(name="DAM Score", 
                          low = "gray30", 
                          high ="red",
                          limits=c(scale_lims), 
                          oob=scales::squish)
  plot(p)
  
}

dev.off()

pdf("old_wt.pdf")
for(df in md_split[metadata$sample[metadata$age!="7 months" & metadata$genotype=="WT"]]) {
  p=ggplot(df)+ geom_point(aes(x=img_x,y=img_y,col=dam_scr))+coord_fixed()+theme_void()+
    ggtitle(paste0(c(df$sample[1], df$genotype[1], df$age[1], df$sex[1]), collapse = "_"))+
    theme(plot.title = element_text(hjust=0.5))+
    scale_color_gradient2(name="DAM Score", 
                          low = "gray30", 
                          high ="red",
                          limits=c(scale_lims), 
                          oob=scales::squish)
  plot(p)
  
}

dev.off()

pdf("old_app.pdf")
for(df in md_split[metadata$sample[metadata$age!="7 months" & metadata$genotype!="WT"]]) {
  p=ggplot(df)+ geom_point(aes(x=img_x,y=img_y,col=dam_scr))+coord_fixed()+theme_void()+
    ggtitle(paste0(c(df$sample[1], df$genotype[1], df$age[1], df$sex[1]), collapse = "_"))+
    theme(plot.title = element_text(hjust=0.5))+
    scale_color_gradient2(name="DAM Score", 
                          low = "gray30", 
                          high ="red",
                          limits=c(scale_lims), 
                          oob=scales::squish)
  plot(p)
}

dev.off()
temp2=md_ser_all %>% group_by(age, genotype, anot) %>% summarise(dam_scr=mean(dam_scr))

temp2$age=gsub("14 months", "old", temp2$age)
temp2$age=gsub("7 months", "young", temp2$age)
temp2$genotype=gsub("APP23", "APP", temp2$genotype)
temp3=temp2 %>% pivot_wider(id_cols="anot", names_from = c( genotype,age), values_from = dam_scr)
temp3$anot=gsub("Olfactory Tubercle","Amygdala",temp3$anot)
md=readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/figures/spat_plots_generic/spat_data_030B.rds")

md=left_join(md, temp3
             )
md$pos1=md$imagecol
md$pos2=md$imagerow
p1=ggplot(md)+ geom_point(aes(x=pos1,y=pos2, col=WT_young),size=0.05)+theme_void()+
  ggtitle(glue("NTG 7 Mo."))+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_fixed() + scale_x_reverse()+scale_y_reverse()+
  scale_color_gradient2(low = "blue", mid = "white", high="red", 
                        limits=c(-2,2),oob=scales::squish)+coord_fixed()


p2=ggplot(md)+ geom_point(aes(x=pos1,y=pos2, col=APP_young),size=0.05)+theme_void()+
  ggtitle(glue("APP23-TG 7 Mo."))+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_fixed() + scale_x_reverse()+scale_y_reverse()+coord_fixed()


p3=ggplot(md)+ geom_point(aes(x=pos1,y=pos2, col=WT_old), size=0.05)+theme_void()+
  ggtitle(glue("NTG 14 Mo."))+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_fixed() + scale_x_reverse()+scale_y_reverse()+coord_fixed()

p4=ggplot(md)+ geom_point(aes(x=pos1,y=pos2, col=APP_old), size=0.05)+theme_void()+
  ggtitle(glue("APP23-TG 14 Mo."))+
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        panel.grid = element_blank())+
  coord_fixed() + scale_x_reverse()+scale_y_reverse()+coord_fixed()

p5=(((p1/p3)|(p2/p4))&scale_color_gradient2(name="Activated Microglia Score", 
                                         low = "gray45" ,
                                         mid ="gray95",
                                         high="red", limits=c(min(temp2$dam_scr),max(temp2$dam_scr)), oob=scales::squish))+
  plot_layout(guides="collect")& theme(legend.position = 'bottom')

for(grp in colnames(temp3)[2:5]) {
  values = temp3[[grp]]
  
  # Define the range and midpoint of your scale
  min_val = min(temp2$dam_scr)
  max_val = max(temp2$dam_scr)
  mid_val = 0
  
  # Adjust the gradient function to correctly map the mid value to 0
  gradient_fn = gradient_n_pal(c("gray45", "gray95", "red"), 
                                values = rescale(c(min_val, mid_val, max_val)))
  
  # Rescale the values to [0, 1] based on your limits
  scaled_values = scales::rescale(values, to = c(0, 1), from = c(min_val, max_val))
  
  # Apply the gradient function to your scaled values
  colors = gradient_fn(scaled_values)
  
  # Show the resulting colors
  print(colors)
  scales::show_col(colors)
  temp4=temp3
  temp4$colr=colors
  system(glue::glue("cp ~/desp1/precast/prec_c25q25g3000/svg_plotting/to_rep.svg ",
                    "~/desp1/precast/prec_c25q25g3000/svg_plotting/dam_scrs/{grp}.svg"))
  i=1
  for(oc in fils2){
    system(glue::glue("sed -i 's/{oc}/{temp4$colr[i]}/g' ~/desp1/precast/prec_c25q25g3000/svg_plotting/dam_scrs/{grp}.svg"))
    i=i+1}
}

