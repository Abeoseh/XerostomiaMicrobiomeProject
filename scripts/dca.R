.libPaths( c( .libPaths(), "~/my_R_libs4.5.1") )

suppressPackageStartupMessages(library(dcurves))
suppressPackageStartupMessages(library(gtsummary))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readxl))

## set working directory
getwd()
setwd("/projects/mougeots_research/abeoseh_flemister/XerostomiaMicrobiomeProject")
getwd()

## set output dir
output="./output/without_men"


## read file 
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx")
microbiome_data <- relocate(microbiome_data, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()

## encode group
microbiome_data <- microbiome_data |> mutate(Group_encoded = case_when(Group == "Xerostomic" ~ 1,
                                                                       Group == "Control" ~ 0) ) |> relocate(Group_encoded)


microbiome_data = select(microbiome_data, -c(Condition, Site, "Other;Other", Taxa, Group))


mod <- glm(Group_encoded ~ ., microbiome_data, family = binomial) %>% broom::augment(newdata = microbiome_data, type.predict = "response")

write.csv(mod, paste(output,"/dca/logistic_regression.csv", sep=""))

png(paste(output, "/dca/results.png",sep=""), width = 1000 )
p = dca(Group_encoded ~ ., microbiome_data) %>% 
  plot(smooth = TRUE) +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )
  # plot(smooth = TRUE, cex.main=10, cex.lab=10, cex=10,
  #      font.lab = 2, font.axis = 2, font.main = 2, font.sub = 2)
print(p)
dev.off()


