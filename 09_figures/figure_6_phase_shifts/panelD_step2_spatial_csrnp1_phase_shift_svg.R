library(tidyverse)
library(patchwork)
library(scico)
library(scales)


setwd("/home/agelber/desp1/precast/prec_c25q25g3000/svg_plotting/csrnp1_phase")

csrnp1_values=readRDS("csrnp1_vals_new.rds") %>% dplyr::filter(age!="young")


csrnp1_values$cluster=factor(csrnp1_values$cluster, 
                             levels = names(readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/svg_plotting/colors.rds")))
csrnp1_values =csrnp1_values%>% arrange(cluster)

mm <- 12
k  <- 0.80  
brks <- c(-mm, -mm*k, 0,  mm*k,  mm)
cols <- c("#000000",
          grDevices::hcl(h = 230, c = 70, l = 45),  # deep blue
          "#FFFFFF",
          grDevices::hcl(h =  10, c = 70, l = 45),  # maroon
          "#000000")

# Create a function to map the colors
gradient_fn <- colorRampPalette(cols)
# Define the rangg
min_val <- -12
max_val <- 12

# Rescale the values to [0, 1] based on your limits
scaled_values <- scales::rescale(csrnp1_values$phi_dif, to = c(0, 1), from = c(min_val, max_val))

# Apply the gradient function to your scaled values
colors=letters[1:23]
for( i in 1:length(colors)) {
  if(scaled_values[i] ==0) {
    colors[i]=gradient_fn(100000)[1]
  } else {
    colors[i]= gradient_fn(100000)[round(scaled_values[i]*100000)]
  }
}

# Show the resulting colors
print(colors)
scales::show_col(colors)
system(glue::glue("cp ~/desp1/precast/prec_c25q25g3000/svg_plotting/to_rep.svg ",
                  "~/desp1/precast/prec_c25q25g3000/svg_plotting/csrnp1_phase/old_new_pal.svg"))
i=1
for(oc in readRDS("/home/agelber/desp1/precast/prec_c25q25g3000/svg_plotting/colors.rds")){
  system(glue::glue("sed -i 's/{oc}/{colors[i]}/g' ~/desp1/precast/prec_c25q25g3000/svg_plotting/csrnp1_phase/old_new_pal.svg"))
  i=i+1}
