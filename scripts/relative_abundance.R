suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(legendry))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(writexl))



## Rscript ./scripts/relative_abundance.R "./csv_files/associated_timepoint/lognorm_data.csv" "./output/associated_timepoint/relative_abundance_" "hand associated" floor
## Rscript ./scripts/relative_abundance.R "./csv_files/skinassoc_timepoint/lognorm_data.csv" "./output/skinassoc_timepoint/relative_abundance_" hand "hand associated"
## Rscript ./scripts/relative_abundance.R "./csv_files/skin_floor_timepoint/lognorm_data.csv" "./output/skin_floor_timepoint/relative_abundance_" hand floor




print("started ")
getwd()
setwd("C:/Users/brean/Downloads/masters/Atrium/")
getwd()

#### set output ####
# output="./output/without_men_10percent_filtered/relative_abundance/relative_abundance.png"
# output="./output/without_men_saliva/relative_abundance/relative_abundance.png"
# output="./output/without_men_saliva_10percent_filtered/relative_abundance/relative_abundance.png"


#### read files ####
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva.xlsx")
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx")

microbiome_data$Group[microbiome_data$Group == "Control"] <- "Non-Xerostomia"

# lognorm_data$timepoint <- paste(lognorm_data$timepoint, "-",lognorm_data$Phenotype, "-" ,lognorm_data$sample_name, sep="")
# 
# plot_data <- lognorm_data %>% select(-c(sample_name, Study_ID, Phenotype, case, ID))  %>% 
#   pivot_longer(!timepoint, names_to = "bacteria", values_to = "count") 


plot_data <- microbiome_data[,1:length(microbiome_data)-1] %>% select(-c(`Other;Other`, Condition, Site))  %>%
  pivot_longer(-c(Taxa, Group), names_to = "bacteria", values_to = "count", )




#### get the top 20 bacteria ####
top_12 <- aggregate(plot_data$count, list(plot_data$bacteria), FUN=sum) %>% arrange(desc(x)) ## aggregate renames bacteria to Group.1 and count to x
top_12 <- top_12$Group.1[1:12]
plot_data$bacteria[!plot_data$bacteria %in% top_12 ] <- "Other"
plot_data$bacteria <- gsub(";s__"," ", gsub("g__", "", plot_data$bacteria))
plot_data$bacteria <- relevel(as.factor(plot_data$bacteria), "Other")




# default_colors <- scales::hue_pal()(length(unique(plot_data$bacteria)))  # ggplot default palette


default_colors <- c("#808080", "#ffa07a", "#E18A00", "#0000FF", "cornflowerblue", 
                    "#FF00FF", "#FFA500", "#8000FF", "#e6e6fa", 
                    "#00FF7F", "#FFD700", "forestgreen", "#FF0000")

names(default_colors) <- unique(plot_data$bacteria)

#### plot the data ####


png(output, width = 1600)


plot <- ggplot(plot_data) +
  
  geom_bar(aes(x = Taxa, y = count, fill = bacteria, color = bacteria), stat = "identity", position = "fill") +
  
  labs(title = "Relative Abundance", x = "Sample Type", y = "") +
  
  facet_grid_paginate(. ~ Group, scales = "free") +
  scale_color_manual(values = default_colors) +
  scale_fill_manual(values = default_colors) +
  
  theme_grey( base_size = 18 ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom",
        strip.text = element_text(face = "bold")) 


print(plot)
dev.off()
  


print("script complete")
