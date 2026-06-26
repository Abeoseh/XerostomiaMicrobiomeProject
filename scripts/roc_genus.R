.libPaths( c( .libPaths(), "~/my_R_libs") )

library(pROC)
library(combiroc)
library(readxl)
library(dplyr)
## set working directory
getwd()
#setwd("X:/OralMedRsch/Abeoseh Flemister/XerostomiaMicrobiomeProject/")
#getwd()

# setwd("C:/Users/brean/Downloads/masters/Atrium/Xerostomia")


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



args = commandArgs(trailingOnly = TRUE)
microbiome_data = read_excel(args[1])
output = args[2]


microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()

metadata <- microbiome_data[,1:4]
metadata$Condition[metadata$Condition == "None"] <- ""
metadata$Group_Condition <- paste(metadata$Group, metadata$Condition)
metadata$Group_Condition_Site <- paste(metadata$Group_Condition, "-", metadata$Site, sep = "")
metadata$Group_Site = paste(metadata$Group, "-", metadata$Site, sep = "")

data_cols <- microbiome_data[,7:length(microbiome_data)-1] %>% t %>% as.data.frame()
data_cols$genus <- gsub("g__", "", gsub(";s__.*", "", rownames(data_cols)))
data_cols <- aggregate(. ~ genus, data = data_cols, FUN=sum) %>% t()

colnames(data_cols) <- data_cols[1,]
microbiome_data <- cbind(metadata, data_cols[2:nrow(data_cols),]) 


microbiome_data <- microbiome_data |> mutate(Group_encoded = case_when(Group == "Xerostomic" ~ 1,
                                                    Group == "Control" ~ 0) ) |> relocate(Group_encoded)
above_0.7_df = data.frame()
auc.df = data.frame()
for (i in 9:length(microbiome_data)){
  bacteria = colnames(microbiome_data)[i]
  bacteria_data = microbiome_data[[bacteria]] %>% as.numeric
  # print(bacteria)
  # print(bacteria_data)
  
  roc_curve = roc(response = microbiome_data$Group_encoded, predictor = bacteria_data)
  
  thres_spec_sens = data.frame("Thres." = roc_curve$thresholds, "Spec." = roc_curve$sensitivities, "Sens." = roc_curve$specificities)
  thres_spec_sens <- thres_spec_sens[thres_spec_sens$Thres != -Inf,]
  thres_spec_sens <- filter(thres_spec_sens, Thres. == min(thres_spec_sens$Thres))
  
  current_auc_df = data.frame(Taxa = bacteria, thres_spec_sens, AUC = auc(roc_curve))  
  auc.df = rbind(auc.df, current_auc_df)
  
  if(auc(roc_curve) > 0.7 | bacteria %in% c("Haemophilus", "Porphyromonas", "Streptococcus")){
    # print("entered confusion matrix")
    # print(bacteria)
    # 2. Extract matrix counts at the best threshold
    counts <- coords(roc_curve, x = thres_spec_sens$Thres, ret = c("threshold", "tp", "tn", "fp", "fn"))[1,]
    # 3. Format into a standard 2x2 matrix

    confusion_matrix <- matrix(
      c(counts$tp, counts$fp, counts$fn, counts$tn), 
      nrow = 2, 
      byrow = TRUE,
      dimnames = list(Predicted = c("Predicted Xerostomia", "Predicted Non-Xerostomia"), Actual = c("Actual Xerostomia", "Actual Non-Xerostomia"))
      ) 
    
    current_0.7_df <- data.frame(Taxa = bacteria, thres_spec_sens, AUC = auc(roc_curve),
                                 "Correctly Predicted Xerostomia" = counts$tp/filter(count(microbiome_data, Group), Group == "Xerostomic")$n, 
                                 "Correctly Predicted Non-Xerostomia" = counts$tn/filter(count(microbiome_data, Group), Group == "Control")$n, 
                                 check.names = FALSE)  
    above_0.7_df = rbind(above_0.7_df, current_0.7_df)
    write.csv(confusion_matrix, paste(output,"/roc_genus/confusion_matrix/", bacteria, ".csv", sep = ""))
    
  }
    
  
  ## remove * from the bacteria names
  bacteria = gsub("\\*", "", bacteria)
  
  png(paste(output, "/roc_genus/curves/ROC for ",bacteria,".png", sep=""))
  plot(roc_curve, main = paste("ROC for", bacteria), print.auc=TRUE, col = "blue")
  dev.off()
  
  
}





write.csv(auc.df, paste(output,"/roc_genus/auc_df.csv",sep=""), row.names = FALSE)


#### CombiROC on Streptococcus, Phrophymonas, and Hemophillus since those were the predictive bacteria
## https://cran.r-project.org/web/packages/combiroc/vignettes/combiroc_vignette_1.html
dat_sub <- microbiome_data %>% 
  select(Taxa, Group, Streptococcus, Porphyromonas, Haemophilus) %>%
  mutate(across(c(Streptococcus, Porphyromonas, Haemophilus), as.numeric))

colnames(dat_sub)[2] = "Class" ## response must be called class

## get the best threshold
data_long <- combiroc_long(dat_sub)
distr <- markers_distribution(data_long, case_class = 'Xerostomic', 
                              y_lim = 0.0015, x_lim = 3000, 
                              signalthr_prediction = TRUE, 
                              min_SE = 0, min_SP = 0, 
                              boxplot_lim = 2000)
max = filter(distr$Coord, is.finite(threshold)) %>% slice_max(Youden)


# Run CombiROC

# combithr, instead, should be set exclusively depending on the needed stringency: 1 is the less stringent and most common choice (meaning that at least one marker in a combination needs to reach the threshold
cb <- combi(
  data = dat_sub,
  signalthr = max$threshold, 
  combithr = 1, 
  case_class="Xerostomic",
  max_length = 3
)

needed_combo = which(cb$Markers == "Haemophilus-Porphyromonas-Streptococcus") - 3
reports <- roc_reports(dat_sub, markers_table = cb,
                        case_class = "Xerostomic",
                        selected_combinations = needed_combo)

reports$Metrics$Thres <- max$threshold
reports$Metrics$Taxa <- cb$Markers[row.names(cb) == row.names(reports$Metrics)]

confusion_matrix <- matrix(
  c(reports$Metrics$TP, reports$Metrics$FP, reports$Metrics$FN, reports$Metrics$TN), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(Predicted = c("Predicted Xerostomia", "Predicted Non-Xerostomia"), Actual = c("Actual Xerostomia", "Actual Non-Xerostomia"))
) 

current_0.7_df <- data.frame(Taxa = "Haemophilus, Porphyromonas, & Streptococcus", 
                             "Thres." = max$threshold, "Spec." = reports$Metrics$SP, "Sens." = reports$Metrics$SE, 
                             AUC = auc(roc_curve), 
                             "Correctly Predicted Xerostomia" = reports$Metrics$TP/filter(count(microbiome_data, Group), Group == "Xerostomic")$n, 
                             "Correctly Predicted Non-Xerostomia" = reports$Metrics$TN/filter(count(microbiome_data, Group), Group == "Control")$n, 
                             check.names = FALSE)  
above_0.7_df = rbind(above_0.7_df, current_0.7_df)

write.csv(above_0.7_df, paste(output,"/roc_genus/AUC_confusion_above0.7.csv", sep = ""), row.names = FALSE)
write.csv(reports$Metrics, paste(output,"/roc_genus/AUC_Haemophilus-Porphyromonas-Streptococcus.csv", sep = ""))
write.csv(confusion_matrix, paste(output,"/roc_genus/confusion_matrix/Haemophilus-Porphyromonas-Streptococcus.csv", sep = ""))
