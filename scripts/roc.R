.libPaths( c( .libPaths(), "~/my_R_libs") )

library(pROC)
library(combiroc)
library(readxl)
library(dplyr)
## set working directory
getwd()
#setwd("X:/OralMedRsch/Abeoseh Flemister/XerostomiaMicrobiomeProject/")
#getwd()


## set output dir
# output="./output/without_men"
# output="./output/with_men"
# output="./output/without_men_saliva"
# output="./output/without_men_10percent_filtered"
output="./output/without_men_saliva_10percent_filtered"


## read file  
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx")
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx")


microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()

metadata <- microbiome_data[,1:4]
metadata$Condition[metadata$Condition == "None"] <- ""
metadata$Group_Condition <- paste(metadata$Group, metadata$Condition)
metadata$Group_Condition_Site <- paste(metadata$Group_Condition, "-", metadata$Site, sep = "")
metadata$Group_Site = paste(metadata$Group, "-", metadata$Site, sep = "")



microbiome_data <- microbiome_data |> mutate(Group_encoded = case_when(Group == "Xerostomic" ~ 1,
                                                    Group == "Control" ~ 0) ) |> relocate(Group_encoded)
print(colnames(microbiome_data)[7:10])

auc.df = data.frame()
for (i in 7:length(microbiome_data)){
  bacteria = colnames(microbiome_data)[i]
  bacteria_data = microbiome_data[[bacteria]]
  
  roc_curve = roc(response = microbiome_data$Group_encoded, predictor = bacteria_data)
  
  thres_spec_sens = data.frame("Thres" = roc_curve$thresholds, "Spec." = roc_curve$sensitivities, "Sens." = roc_curve$specificities)
  thres_spec_sens <- thres_spec_sens[thres_spec_sens$Thres != -Inf,]
  thres_spec_sens <- filter(thres_spec_sens, Thres == min(thres_spec_sens$Thres))
  
  current_auc_df = data.frame(Taxa = bacteria, thres_spec_sens, AUC = auc(roc_curve))  
  auc.df = rbind(auc.df, current_auc_df)
  
  ## remove * from the bacteria names
  bacteria = gsub("\\*", "", bacteria)
  
  png(paste(output, "/roc/ROC for ",bacteria,".png", sep=""))
  plot(roc_curve, main = paste("ROC for", bacteria), print.auc=TRUE, col = "blue")
  dev.off()
  
  
}





write.csv(auc.df, paste(output,"/roc/auc_df.csv",sep=""), row.names = FALSE)


# auc_above_70 = filter(auc.df, AUC > 0.70)
# roc_above_70 = select(microbiome_data, Taxa, Group) |> 
#   cbind(microbiome_data[,c(colnames(microbiome_data)[ colnames(microbiome_data) %in% auc_above_70$bacteria ]) ])
# colnames(roc_above_70)[colnames(roc_above_70) == "Group"] = "Class" ## the response must be called Class
# colnames(roc_above_70)[3:length(roc_above_70)] = gsub(";|-", "_", colnames(roc_above_70)[3:length(roc_above_70)]) ## semi-colons not allowed
# colnames(roc_above_70)[3:length(roc_above_70)] = gsub("*", "", colnames(roc_above_70)[3:length(roc_above_70)]) ## semi-colons not allowed
# 
# roc_long <- combiroc_long(roc_above_70)
# 
# distr <- markers_distribution(roc_long, case_class = "Xerostomic",
#                               signalthr_prediction = TRUE,
#                               min_SE = 0, min_SP = 0,
#                               boxplot_lim = max(roc_long$Values) + 0.5)
# distr$Density_plot
# 
# tab = combi(roc_above_70, signalthr = 0, combithr = 1, case_class="Xerostomic", max_length = dim(auc_above_70)[1] ) ## signalthr is from distr$Density_plot
# 
# reports <- roc_reports(roc_above_70, markers_table = tab,
#                        case_class = "Xerostomic",
#                        single_markers = c(row.names(tab)[ 1:dim(auc_above_70)[1] ]),
#                        selected_combinations = c( dim(tab)[1] - dim(auc_above_70)[1] ))
# 
# #print("start report")
# #bacteria_for_combiroc = row.names(tab)[ 1:dim(auc_above_70)[1] ]
# 
# 
# #print(paste((length(bacteria_for_combiroc) + 1), dim(tab)[1] ))
# #print(dim(tab)[1] - dim(auc_above_70)[1])
# 
# reports_csv <- roc_reports(roc_above_70, markers_table = tab,
#                        case_class = "Xerostomic",
#                        single_markers = c(row.names(tab)[ 1:dim(auc_above_70)[1] ]),
#                        selected_combinations = c( 1:(dim(tab)[1] - dim(auc_above_70)[1])  ))
# 
# 
# #print("end report")
# 
# AUC_report <- merge(select(tab, Markers), reports_csv$Metrics, by=0)
# AUC_report$Row.names = NULL
# row.names(AUC_report) = AUC_report$X
# AUC_report$X = NULL
# 
# png(paste(output,"/combiroc/combiroc_plot.png",sep=""))
# reports$Plot
# dev.off()
# 
# write.csv(tab, paste(output,"/combiroc/AUCs.csv",sep=""))
# write.csv(AUC_report, paste(output,"/combiroc/report.csv",sep=""))
# 
# 
# 
