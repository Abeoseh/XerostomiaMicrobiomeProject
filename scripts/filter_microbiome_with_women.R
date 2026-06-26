.libPaths( c( .libPaths(), "~/my_R_libs4.5.1") )

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))

df <- read_excel("../csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
df2 = read.csv("../output/without_men/roc/auc_df.csv")



metadata_cols = select(df, colnames(df)[1:5])
num_col = select(df, "...491")
bacteria = filter(df2, auc >= 0.7)$bacteria
bacteria_df = select(df, all_of(bacteria))
dim(bacteria_df)
colnames(bacteria_df)
bacteria
df3 = cbind(metadata_cols,bacteria_df, num_col)
head(df3)
colnames(df)[1:5]
metadata_cols = select(df, colnames(df)[1:5])
df3 = cbind(metadata_cols,bacteria_df, num_col)
colnames(df3)

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", df3, rowNames = FALSE)
saveWorkbook(wb, "../csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx", overwrite = TRUE)
