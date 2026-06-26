rm(list=ls())
.libPaths( c( .libPaths(), "~/my_R_libs") )
library(ggplot2)
library(caret)
library(pROC)
library(dplyr)
library(rfPermute)
suppressPackageStartupMessages(library(readxl))

args = commandArgs(trailingOnly = TRUE)
microbiome_data = read_excel(args[1])
output = args[2]

#####################################################
setwd("C:/Users/brean/Downloads/masters/Atrium")
microbiome_data = read_excel("Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx")
output = "./output"

#### 10 fold cross-validation #########################
# Load your data
microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()


microbiome_data <- microbiome_data |> mutate(Group_encoded = case_when(Group == "Xerostomic" ~ 1,
                                                                       Group == "Control" ~ 0) ) |> relocate(Group_encoded)
# Define 10-fold cross-validation
microbiome_data$SubjectID <- gsub("\\.(.*)", "", microbiome_data$Taxa)

set.seed(2)
folds <- groupKFold(as.factor(microbiome_data$SubjectID) , k = 10) ## keep all 3 sites together

microbiome_data = select(microbiome_data, -c(Condition, "Other;Other", Group, Taxa))
colnames(microbiome_data) <- make.names(colnames(microbiome_data)) ## remove ; and *
microbiome_data$Group_encoded <- as.factor(microbiome_data$Group_encoded)
microbiome_data <- relocate(microbiome_data, SubjectID)

# valid_folds <- list()
rf_models=list()
rf_predictions=list()


for (i in seq_along(folds)) {
  train_idx <- folds[[i]]
  train_data <- microbiome_data[train_idx, ]
  test_data <- microbiome_data[-train_idx, ]
  # print(train_idx)
  
  # Remove Site from features
  train_data_rf <- train_data %>% select(-SubjectID, -Site)
  test_data_rf <- test_data %>% select(-SubjectID, -Site)
  print(unique(test_data_rf$Group_encoded))
  
  
  model <- rfPermute(Group_encoded ~ ., data = train_data_rf, importance = TRUE, num.rep=500)
  pred <- predict(model, newdata = test_data_rf, type = "prob")

  rf_models[[length(rf_models) + 1]] <- model
  rf_predictions[[length(rf_predictions) + 1]] <- pred
  valid_folds[[length(valid_folds) + 1]] <- test_idx
}

#Get importance and p-values only from valid models
importance_long <- lapply(rf_models, function(model) {
    imp <- importance(model)
    data.frame(
      Taxon = rownames(imp),
      Importance = imp[, "MeanDecreaseGini"],
      pvals=imp[,'MeanDecreaseGini.pval'],
      stringsAsFactors = FALSE )
})

## Make the importance into a df
importance_all <- do.call(rbind, lapply(seq_along(importance_long), function(i) {
  df <- importance_long[[i]]
  df$Fold = paste("Fold:",i) 
  rbind(df, source = rep(names(importance_long)[i], nrow(df)))
}))

importance_all$pvalue_corrected = p.adjust(importance_all$pvals, method = "BH")
write.csv(importance_all, paste(output, "/importance_bars.csv",sep="" ), row.names = F)

importance_all <- importance_all %>%
  mutate(Fold = factor(Fold, levels = unique(Fold[order(as.numeric(gsub("Fold: ", "", Fold)))])) )

g <- ggplot(importance_all, aes(x = reorder(Taxon, -Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = ifelse(importance_all$pvalue_corrected < 0.49,'blue','black')) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  facet_wrap(~Fold) +
  labs(title = "Important Taxa",
       x = "Taxa", y = "Gini Importance")

png(paste(output,"/plots/importance_bars.png", sep=""), width = 600, height=500)
print(g)
dev.off()

roc_df = data.frame()
auc_df = data.frame()
for (i in 1:length(rf_predictions)){
  train_idx <- folds[[i]]
  test_data <- microbiome_data[-train_idx, ] %>% select(-SubjectID, -Site)
  if (length(unique(test_data[,1])) == 2){

    rf_roc <- roc(test_data[,1], rf_predictions[[i]][,1])
    
    roc_df <- rbind(roc_df, data.frame(Fold = paste("Fold:",i), sensitivity = rf_roc$sensitivities, specificity = rf_roc$specificities ))
    auc_df <- rbind(auc_df, data.frame(Fold = i, AUC = auc(rf_roc) ))
  }
}

write.csv(roc_df, paste(output, "/roc_df.csv", sep=""), row.names = F)
write.csv(auc_df, paste(output, "/auc_df.csv", sep=""), row.names = F)


g <- ggplot(auc_df, aes(x = as.factor(Fold), y = AUC)) +
  geom_bar(stat = "identity", fill = "tomato") +
  ylim(0, 1) +
  theme_minimal(base_size = 14) +
  labs(title = "AUC per Fold", x = "Fold", y = "AUC")

png(paste(output, "/plots/AUC_per_fold.png", sep=""))
print(g)
dev.off()







