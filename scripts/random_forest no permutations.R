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



microbiome_data = read_excel(args[1])
output = args[2]


## read file 
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx")



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

## 10 fold CV
# split data into 10 chunks
# train on 9 of the chunks and test on the 10th. Repeat 10 times


## make the partitions
set.seed(100)
k = 10
folds <- createFolds(y = microbiome_data$Group_encoded, k = k, list = TRUE, returnTrain = FALSE)

# print(partitions)
auc_df = data.frame()
aucs = c()
for (i in seq_along(folds)){

  testing = microbiome_data[ folds[[i]], ]

  training = microbiome_data[ setdiff(seq_len(nrow(microbiome_data)), folds[[i]]), ]
  
  set.seed(100)
  RF_fit <- randomForest(Group_encoded~., data = training, importance = TRUE)
  
  set.seed(100)
  RF_pred <- predict(RF_fit, testing, type = "prob")
  rf_roc <- roc(testing[,1], RF_pred[,1])
  
  auc_df = rbind(auc_df, data.frame(CV = i, AUC = auc(rf_roc)))
  aucs = append(aucs, auc(rf_roc))

  ## make the roc curves
  png(paste(output, "/plots/ROC_cv", i, ".png", sep=""))
  p <- plot(rf_roc, col = "red", print.auc = TRUE)
  title(paste("ROC Curve of CV: ", i, sep=""), line = + 2.5, cex.main=1.5)
  p
  dev.off()
  
  
  ## make the feature importance plot
  feat_imp_df <- importance(RF_fit) %>%
    data.frame() %>%
    mutate(feature = row.names(.))
  

  feat_imp_df$feature = row.names(feat_imp_df)
  

  feat_imp_df$abs_MeanDecreaseGini <- abs(feat_imp_df$MeanDecreaseGini)

  feat_imp_df = arrange(feat_imp_df, desc(abs_MeanDecreaseGini))


  feat_imp_top_50 = feat_imp_df[1:min(nrow(feat_imp_df), 25),]


  # plot dataframe
  g <-  ggplot(feat_imp_top_50, aes(x = reorder(feature, abs_MeanDecreaseGini),y = abs_MeanDecreaseGini)) +
    geom_bar(stat='identity') +
    coord_flip() +
    theme_classic() +
    labs(
      x = "Feature",
      y = "Importance",
      title = paste("Variable Importance Plot for CV ", i, sep="") ) +
    theme(plot.title = element_text(size=15), axis.text=element_text(size=11),
          axis.title=element_text(size=15), axis.text.y = element_text(face = "bold"))

  png(paste(output,"/plots/importance_bars_cv", i, ".png", sep=""), width = 600, height=500)
  print(g)
  dev.off()
  
  # do permutations and plot data
  group.rp <- rfPermute(Group_encoded ~ ., data = training, na.action = na.omit, ntree = 100, num.rep = 50)
  print(attributes(group.rp))
  print(group.rp$pval)
  print("-----------------------------------")
  
  g <- plotImportance(group.rp, scale = TRUE)
  
  png(paste(output,"/plots/importance_bars_permutations_cv", i, ".png", sep=""), width = 600, height=500)
  print(g)
  dev.off()
  
  
}

print("AUCs:")
print(aucs)
print(paste("average AUC:", mean(aucs)))

auc_df = rbind(auc_df, data.frame(CV = "avg", AUC = mean(aucs)))
write.csv(auc_df, paste(output,"/AUCs.csv", sep=""), row.names = F)


print("__________________________________________________________________________")
cv = trainControl(method = "cv", number = 10)
RF_fit <- train(Group_encoded~., data = microbiome_data, method = "rf",
                trControl = cv)

print(RF_fit)





