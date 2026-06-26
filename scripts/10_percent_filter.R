## library packages  ####
.libPaths( c( .libPaths(), "~/my_R_libs") )

# suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
# suppressPackageStartupMessages(library(vegan))
# suppressPackageStartupMessages(library(ecodist))
# suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(writexl))
# suppressPackageStartupMessages(library())


##### set working directory ####
# getwd()
# setwd("X:/OralMedRsch/Abeoseh Flemister/XerostomiaMicrobiomeProject/")
# getwd()



#### read file ####
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva.xlsx")
microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)] %>% as.data.frame()

data <- microbiome_data[,7:length(microbiome_data)-1]

percent_nonzero <- colSums(data != 0)/dim(data)[1]
ten_percent_data <- data[,names(percent_nonzero[percent_nonzero > .1])]

ten_percent_df <- cbind(microbiome_data[,1:6], ten_percent_data, microbiome_data["...491"])
colnames(ten_percent_df)[length(ten_percent_df)] <- "...491"

write_xlsx(ten_percent_df, "csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx", col_names = TRUE)

