.libPaths( c( .libPaths(), "~/my_R_libs") )
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(reprtree))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rfPermute))


args = commandArgs(trailingOnly = TRUE)

## set output dir
# output="./output/without_men/RandomForest_ROCabove0.7"

# output="./output/without_men/RandomForest"
# output="./output/with_men/RandomForest"

# Rscript ./scripts/random_forest.R "./csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx" ./output/without_men/RandomForest_ROCabove0.7

microbiome_data = read_excel(args[1])
output = args[2]


## read file 
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
# microbiome_data <- read_excel("Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx")
# output = "./output/my_rf"


microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()

# metadata <- microbiome_data[,1:4]
# metadata$Condition[metadata$Condition == "None"] <- ""
# metadata$Group_Condition <- paste(metadata$Group, metadata$Condition)
# metadata$Group_Condition_Site <- paste(metadata$Group_Condition, "-", metadata$Site, sep = "")
# metadata$Group_Site = paste(metadata$Group, "-", metadata$Site, sep = "")



microbiome_data <- microbiome_data |> mutate(Group_encoded = case_when(Group == "Xerostomic" ~ 1,
                                                                       Group == "Control" ~ 0) ) |> relocate(Group_encoded)

microbiome_data$ID = sapply(strsplit(microbiome_data$Taxa,"\\."), `[`, 1)
microbiome_data = select(microbiome_data, -c(Condition, Site, "Other;Other", Group, Taxa))
microbiome_data = aggregate(. ~ ID, microbiome_data, FUN = sum)
microbiome_data <- microbiome_data |> mutate(Group_encoded = case_when(Group_encoded > 0 ~ 1,
                                                                       Group_encoded == 0 ~ 0) )
microbiome_data$ID <- NULL

# print(dim(microbiome_data))
# print(head(microbiome_data))

colnames(microbiome_data) <- make.names(colnames(microbiome_data)) ## remove ; and *


microbiome_data$Group_encoded <- as.factor(microbiome_data$Group_encoded)


# colnames(microbiome_data)
# colnames(microbiome_data)[1:10]

### 10 fold CV
## split data into 10 chunks
## train on 9 of the chunks and test on the 10th. Repeat 10 times


## make the partitions
set.seed(100)
k = 10
folds <- createFolds(y = microbiome_data$Group_encoded, k = k, list = TRUE, returnTrain = FALSE)

# print(partitions)
auc_df = data.frame()
roc_df = data.frame()
aucs = c()
rf_models=list()

for (i in seq_along(folds)){

  testing = microbiome_data[ folds[[i]], ]

  training = microbiome_data[ setdiff(seq_len(nrow(microbiome_data)), folds[[i]]), ]
  
  
  # RF_fit <- randomForest(Group_encoded~., data = training, importance = TRUE)
  
  set.seed(100)
  RF_fit <- rfPermute(Group_encoded ~ ., data = training, importance = TRUE, num.rep=500)
  rf_models[[length(rf_models) + 1]] <- RF_fit
  
  ## predictions
  set.seed(100)
  pred <- predict(RF_fit, newdata = testing, type = "prob")
  
  ## compute ROC and AUC
  rf_roc <- roc(testing[,1], pred[,1])
  
  auc_df = rbind(auc_df, data.frame(CV = i, AUC = auc(rf_roc)))
  aucs = append(aucs, auc(rf_roc))

  roc_df <- rbind(roc_df, data.frame(CV = paste("Fold:",i), specificity = rf_roc$specificities, sensitivity = rf_roc$sensitivities ))
  

  ## plot the roc curves
  png(paste(output, "/plots/ROC_cv", i, ".png", sep=""))
  p <- plot(rf_roc, col = "red", print.auc = TRUE)
  title(paste("ROC Curve of CV: ", i, sep=""), line = + 2.5, cex.main=1.5)
  p
  dev.off()
  
  
}

## make the AUC df and plot the AUCs
print("AUCs:")
print(aucs)
print(paste("average AUC:", mean(aucs)))

auc_df = rbind(auc_df, data.frame(CV = "avg", AUC = mean(aucs)))
write.csv(auc_df, paste(output,"/AUCs.csv", sep=""), row.names = F)
write.csv(roc_df, paste(output,"/ROCs.csv", sep=""), row.names = F)

num_levels <- sort(as.numeric(auc_df$CV[auc_df$CV != "avg"]))
levels <- c(as.character(num_levels), "avg")
auc_df$CV <- factor(auc_df$CV, levels = levels)

avg_auc_line = auc_df$AUC[auc_df$CV == "avg"]

g <- ggplot(filter(auc_df, CV != "avg"), aes(x = CV, y = AUC)) + 
  geom_bar(stat = "identity", fill = "tomato") +
  ylim(0, 1) +
  geom_abline(slope = 0, intercept = avg_auc_line, col = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "AUC per Fold", x = "Fold", y = "AUC")

png(paste(output, "/plots/AUC_per_fold.png", sep=""))
print(g)
dev.off()



### plot the variable importances
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

png(paste(output,"/plots/importance_bars.png", sep=""), width = 600, height=800)
print(g)
dev.off()







