########################################################################################
## From the bacteria which were found to have a AUC > 0.7 with ROC plot the abundance ##
########################################################################################
.libPaths( c( .libPaths(), "~/my_R_libs4.5.1") )

library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(forcats)


args = commandArgs(trailingOnly = TRUE)

input <- args[1]
auc_path <- args[2]
output_dir <- args[3]

# output_dir = "./output/without_men_saliva_10percent_filtered/ROC_fisher"
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx")
# auc_df <- read.csv("./output/without_men_saliva_10percent_filtered/roc/auc_df.csv") %>% filter(AUC > 0.7)



output_dir = "./output/without_men_10percent_filtered/ROC_fisher"
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx")
auc_df <- read.csv("./output/without_men_10percent_filtered/roc/auc_df.csv") %>% filter(AUC > 0.7)



microbiome_data <- read_excel(input)
auc_df <- read.csv(auc_path) %>% filter(AUC > 0.7)


df = select(microbiome_data, Group, all_of(auc_df$Taxa)) %>%
  pivot_longer(auc_df$Taxa, names_to = "Taxa", values_to = "Amount")

df$Group[df$Group == "Control"] <- "NX"
df$Group[df$Group == "Xerostomic"] <- "XS"
df$Group <- factor(df$Group, levels = c("XS", "NX"))

df$Taxa <- gsub("g__|f__|o__|c__|p__|k__|NA ", "", gsub(";s__|;", " ", df$Taxa))
df <- df %>% mutate(Taxa = fct_reorder(Taxa, Amount, .fun = max, .desc = TRUE))


meandf <- df %>%
  # group_by(Taxa, Group) %>%
  group_by(Taxa) %>%
  summarise(means = mean(Amount), .groups = "drop")

###### fisher test: testing association between ##
## Group and amount of samples below and above the mean ######

AboveBelow <- df %>%
  left_join(meandf, by = "Taxa") %>%
  # mutate(Group = coalesce(Group.x, Group.y)) %>%
  # select(-Group.x, -Group.y) %>%
  mutate(
    AboveMean = case_when(
      Amount >= means ~ TRUE,
      Amount < means ~ FALSE,
    )
    ) %>%
  # group_by(Taxa) %>% count(AboveMean)
  group_by(Group,Taxa) %>% count(AboveMean)


results <- AboveBelow %>%
  group_by(Taxa) %>%
  group_modify(~{
    tab <- xtabs(n ~ Group + AboveMean, data = .x)
    ft <- fisher.test(tab)
    
    tibble(
      p.value = ft$p.value,
      odds.ratio = unname(ft$estimate)
    )
  })


results$FDR <-  p.adjust(results$p.value, method = "BH")


write.csv(results, paste(output_dir, "/fisher_results_one_mean.csv",sep=""), row.names = FALSE)


###### plot ###### 


pvals <- results %>%
  mutate(
    label = paste0("Fisher p = ", signif(p.value, 3))
  )
pvals <- pvals %>%
  left_join(
    df %>%
      group_by(Taxa) %>%
      summarise(y = max(Amount) * 0.95),
    by = "Taxa"
  ) %>%
  mutate(x = 0.5) 


p <- ggplot(df, aes(as.factor(Group), Amount)) +
  geom_point() +
  geom_text(data = results,
    aes(x = 1.5, y = 0, label = paste("p-value",sprintf("%.2e", FDR))), inherit.aes = FALSE, size = 3) +
  geom_hline(data = meandf, aes(yintercept = means), color = "cornflowerblue", linetype = "solid", inherit.aes = FALSE) +
  # geom_crossbar(data = meandf,aes(x = as.factor(Group), y = means, ymin = means, ymax = means),
  #   width = 0.6, color = "cornflowerblue", fatten = 0, inherit.aes = FALSE) +
  labs(x = "") +
  facet_wrap(~ Taxa, scales = "free_y")


png(paste(output_dir,"/Bacteria_one_mean.png", sep=""), width = 1000, height = 1000)
print(p)
dev.off()








