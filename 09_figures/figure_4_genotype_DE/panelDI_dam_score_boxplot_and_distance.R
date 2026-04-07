library(tidyverse)
library(scattermore)
library(scales)
library(lmerTest)
library(lme4)
library(nlme)

setwd("~/desp1/precast/prec_c25q25g3000/figures/fig4_app_basic/dam_scr")
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
md_ser_all=left_join(md_ser_all,coords, by=c("sample","cell"))
md_ser_all=md_ser_all %>% ungroup() %>% 
  mutate(dam_scr=-1*(PC1+PC2))


scale_lims=c(min(md_ser_all$dam_scr),max(md_ser_all$dam_scr))

scale_lims=c(quantile(md_ser_all$dam_scr,.001),quantile(md_ser_all$dam_scr,.999))

data_dir="/cndd2/agelber/hal/qc_aligned"

metadata = data.table::fread("/cndd2/agelber/hal/metadata.csv") %>% 
  mutate(sample=str_split_fixed(sample, "_", 2)[,1]) %>% 
  dplyr::filter(sample %in% list.files(data_dir))  %>% arrange(sex)

short_names=readRDS("../../../objects/short_names.rds")

temp=md_ser_all %>% group_by(age, genotype, sex,anot, sample) %>% 
  summarise(dam_scr=mean(dam_scr))%>% mutate(anot=factor(plyr::mapvalues(anot,
                                                                         short_names,
                                                                         names(short_names)) , 
                                                         levels=names(short_names)),
                                             age=factor(age, levels=c("7 months", "14 months")),
  ) %>% dplyr::filter(grepl("CTX[0-9]", anot))
temp$genotype=gsub("APP23", "APP23-TG", temp$genotype)
temp$genotype=gsub("WT", "NTG", temp$genotype)
temp$anot=gsub("Olfactory Tubercle", "Amygdala", temp$anot)
temp2=temp %>% ungroup()%>%
  group_by(age, genotype, anot) %>% 
  summarise(mean=mean(dam_scr), sd=sd(dam_scr), se=sd/sqrt(n()))


genotype_colors=ggsci::pal_nejm()(2)
names(genotype_colors)=c("APP23-TG", "NTG")


ggplot(temp2 ,
       aes(x=anot, y=mean, fill=genotype,
           ymin = mean - se, ymax = mean + se))+
  geom_col(position=position_dodge(width = 0.6), width = 0.6)+
  geom_errorbar(position=position_dodge(width = 0.6), width=0.3 ) +
  geom_point(data=temp,aes(x=anot, y=dam_scr, group=genotype),
             position=position_dodge(width = 0.6),
             size=0.2,
             inherit.aes = F)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_fill_manual(name="Genotype", values=genotype_colors)+
  scale_color_manual(name="Genotype", values=genotype_colors)+
  xlab("")+ylab("DAM Score")+facet_wrap(~age)+scale_y_continuous(limits = c(-1.5,5), oob=scales::squish)




g1=lapply(split(temp, paste0(temp$anot, "_", temp$genotype)),
          function(x) {
            model = lme(dam_scr~ age, random = ~1|sex,
                        data=x)
            
            av=anova(model)
            av$`p-value`[2] }) %>% unlist() %>% sort()

g1

g2=lapply(split(temp, paste0(temp$anot, "_", temp$age)),
          function(x) {
            model = lme(dam_scr~ genotype, random = ~1|sex,
                        data=x)
            
            av=anova(model)
            av$`p-value`[2] }) %>% unlist() %>% sort()
g2

g= lapply(split(temp, paste0(temp$anot, "_", temp$age)),
          function(x) {
            t.test(x$dam_scr[x$genotype!="NTG"], 
                   x$dam_scr[x$genotype=="NTG"], 
                   alternative="g")[["p.value"]]}) %>% unlist() %>% sort()

p.adjust(g)
g

g2= lapply(split(temp, paste0(temp$anot, "_", temp$genotype)),
           function(x) {
             t.test(x$dam_scr[x$age!="7 months"], 
                    x$dam_scr[x$age=="7 months"], 
                    alternative="g")[["p.value"]]}) %>% unlist() %>% sort()

p.adjust(g2)
g2

