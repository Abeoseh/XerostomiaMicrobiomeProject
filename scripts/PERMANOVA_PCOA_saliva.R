## library packages  ####
.libPaths( c( .libPaths(), "~/my_R_libs") )

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ecodist))
suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library())


##### set working directory ####
# getwd()
# setwd("X:/OralMedRsch/Abeoseh Flemister/XerostomiaMicrobiomeProject/")
# getwd()

#### set output dir ####
# output="./output/without_men_saliva"
output="./output/without_men_saliva_10percent_filtered"


#### read files ####
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva.xlsx")
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx")

microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()

microbiome_data$Group[ grep("Control", microbiome_data$Group) ] <- "NX"
microbiome_data$Group[ grep("Xerostomic", microbiome_data$Group) ] <- "XS"
microbiome_data$Group <- factor(microbiome_data$Group, levels = c("XS", "NX"))

metadata <- microbiome_data[,1:4]
metadata$Condition[metadata$Condition == "None"] <- ""
metadata$Group_Condition <- paste(metadata$Group, metadata$Condition)
metadata$Group_Condition_Site <- paste(metadata$Group_Condition, "-", metadata$Site, sep = "")
metadata$Group_Site = paste(metadata$Group, "-", metadata$Site, sep = "")
data_cols <- microbiome_data[,c(6:length(microbiome_data))]

#### PERMANOVAs ####


## PERMANOVA with effect of Group and Site independently
permanova_independent <- adonis2(data_cols ~ Group, data = metadata, method = "bray", by = "terms") %>% as.data.frame()

## PERMANOVA with interaction between Group-Condition and Site for xerostomia groups separated 
permanova_interaction_group_condition <- adonis2(data_cols ~ Group_Condition, data = metadata, method = "bray", by = "terms") %>% as.data.frame()
permanova_interaction_group_condition$method <- "Group_Condition"

## write data to an excel file
wb <- createWorkbook()

addWorksheet(wb, "Group independently")
writeData(wb, "Group independently", permanova_independent, rowNames = TRUE)


saveWorkbook(wb, paste(output,"/permanova/permanova_results.xlsx",sep=""), overwrite = TRUE)



# write.csv(rbind(permanova_independent, permanova_interaction) ,"./output/permanova/permanova_results.csv")


#### PCOA ####
pcoa <- function(df, metadata, chosen_title, metadata_column, pc = NULL){
  
  # View(df)
  # Sys.wait(30)
  # View(metadata)
  bray <- vegdist(df, method = "bray")
  # print("done")
  pcoa_val <- pco(bray, negvals = "zero", dround = 0)

  if(!is.null(pc)){
    
    cat("example to access elements in the list: \n
    l = pcoa(df)\n
    l[1] # returns bray-curtis values\n
    l[2] # returns the pcoa results\n
    l[2][[1]] # returns vectors and values which can be accessed via $\n
    ex: l[2][[1]]$vectors")
    
    return(list(bray, pcoa_val))
  }
  
  pcoa_val.df = data.frame(Taxa = row.names(pcoa_val$vectors),
                           PCoA1 = pcoa_val$vectors[,1], 
                           PCoA2 = pcoa_val$vectors[,2])
  
  pcoa_val.df = merge(pcoa_val.df, metadata, by.all = "Taxa", all.x = TRUE)
  pcoa_val.df[[metadata_column]] <- as.factor(pcoa_val.df[[metadata_column]])
  
  # print(dim(pcoa_val.df))
  # print(pcoa_val.df)
  # print(unique(pcoa_val.df[[metadata_column]]))
  
  # compute variance explained by each PCoA via eigenvalues 
  eigenvalues = pcoa_val$values ## summing the eigenvalues gives total variance since one eigenvalue tells you how much variance is captured by one PCo

  metadata_column <- sym(metadata_column) ## make a symbol for plotting
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = !!metadata_column), alpha = 0.7) +
    # scale_color_manual(values = group.colors) +
    stat_ellipse(level = 0.95, aes(group = !!metadata_column, color = !!metadata_column)) +
    theme(plot.margin = margin(10, 10, 20, 10), plot.caption = element_text(face = "bold", hjust = 0), text=element_text(size=18)) +
    labs(title = chosen_title, x = paste("PCo1 (", round((eigenvalues[1] / sum(eigenvalues)) * 100, 2), "%)",sep=""), 
         y = paste("PCo2 (", round((eigenvalues[2] / sum(eigenvalues)) * 100, 2),"%)",sep=""),
         caption = paste("PERMANOVA p-value =", signif(perm_df[["Pr(>F)"]],2), "; R2 =", round(perm_df$R2, 2)) )
  return(pco.plot)
}

row.names(data_cols) <- metadata$Taxa

png(paste(output,"/permanova/pcoa_group.png",sep=""))
perm_df <- permanova_independent["Group",]
pcoa(data_cols, metadata, "Group", "Group")
dev.off()


png(paste(output,"/permanova/pcoa_group_condition.png",sep=""))
perm_df <- permanova_interaction_group_condition["Group_Condition",]
pcoa(data_cols, metadata, "Group-Condition", "Group_Condition")
dev.off()

used_groups = c()
pairwise_pcoas = data.frame()
for ( group1 in unique(metadata$Group_Condition) ){
  for ( group2 in unique(metadata$Group_Condition) ){
    if ( group1 != group2 & !(group2 %in% used_groups) ){
      print(paste(group1, group2))
      filtered_metadata <- metadata[metadata$Group_Condition == group1 | metadata$Group_Condition == group2,]
      row.names(filtered_metadata) <- filtered_metadata$Taxa
      filtered_data_cols <- data_cols[ c(row.names(filtered_metadata)), ]
      # data_cols <- microbiome_data[,c(6:length(microbiome_data))]
      # metadata <- microbiome_data[,1:4]
      permanova_group_condition <- adonis2(filtered_data_cols ~ Group_Condition, data = filtered_metadata, method = "bray", by = "terms") %>% as.data.frame()
      permanova_group_condition$group1 = group1
      permanova_group_condition$group2 = group2
      
      png(paste(output,"/permanova/pcoa_group_condition", group1, "-", group2, ".png",sep = ""), width=600, height=600)
      perm_df <- permanova_group_condition["Group_Condition",]
      pcoa(filtered_data_cols, filtered_metadata, paste("Group-Condition for", group1, "and", group2), "Group_Condition") %>% print()
      dev.off()
      
      pairwise_pcoas <- rbind(pairwise_pcoas, permanova_group_condition)
      
      
      used_groups = append(used_groups, group1)
  

    }  
}}

addWorksheet(wb, "pairwise Group_Condition")
writeData(wb, "pairwise Group_Condition", pairwise_pcoas, rowNames = TRUE)


saveWorkbook(wb, paste(output,"/permanova/permanova_results.xlsx",sep=""), overwrite = TRUE)




#### alpha diversity ####
alpha_diversity <- function(df, metadata, metadata_column){
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(ggsignif))
  
  # df <- data_cols

  evenness <- diversity(df, "simpson") / log(specnumber(df)) ## calculates evenness index
  
  shannon_diversity <- diversity(df, index = "shannon") ## shannon diversity
  
  simpson_diversity <- diversity(df, "simpson")  ## calculates simpson diversity
  
  alpha_diversity = cbind(metadata, shannon_diversity, evenness, simpson_diversity)
  
  # simpson_diversity |> print()
  # break
  
  rm(evenness, shannon_diversity, simpson_diversity)
  alpha_diversity_tidy <- alpha_diversity |> 
    # select(-c(Group_Condition_Site, Group_Site)) |>
    mutate(Taxa = row.names(df)) |>
    gather(key = alphadiv_index, value = alphadiv_values, 
           -Taxa, -Group, -Condition, -Site, -Group_Condition, -Group_Condition_Site, -Group_Site)
  
  
  # return(alpha_diversity_tidy)
  

  ## plot indices 
  alpha_indices <- unique(alpha_diversity_tidy$alphadiv_index)
  groups <- unique(alpha_diversity_tidy[[metadata_column]])
  
  metadata_column_var <- sym(metadata_column)
  for (index in alpha_indices){
    used_groups = c()
    for (group1 in groups){
      for (group2 in groups){
        if (group1 != group2 & !(group2 %in% used_groups) ){
        used_groups = append(used_groups, group1)
        ## calculate p-values, make the p-values into a df, and only keep significant p-values
        pval <- wilcox.test(alpha_diversity_tidy$alphadiv_values[ alpha_diversity_tidy$alphadiv_index  == index & alpha_diversity_tidy[[metadata_column]] == group1],
                            alpha_diversity_tidy$alphadiv_values[ alpha_diversity_tidy$alphadiv_index  == index & alpha_diversity_tidy[[metadata_column]] == group2])$p.value %>%
                signif(3)
    

        # current_index <- filter(alpha_diversity_tidy, alphadiv_index == !!index & (metadata_column == !!group1 | metadata_column == !!group2 ))
        current_index <- alpha_diversity_tidy[alpha_diversity_tidy$alphadiv_index == index & (alpha_diversity_tidy[[metadata_column]] == group1 | alpha_diversity_tidy[[metadata_column]] == group2 ),]
        
        current_index$Group <- factor(current_index$Group, levels = c("XS", "NX"))
        alpha_plot = ggplot(current_index, aes(!!metadata_column_var, alphadiv_values)) +
          labs(title = gsub("_", " ", index), y = index, caption = paste("p-value:", pval)) +
          geom_boxplot(aes(color = !!metadata_column_var)) +
          geom_point(aes(color = !!metadata_column_var)) +
          theme(axis.text.x = element_text(face = "bold", vjust = 0.5, hjust=1, size=30), legend.position="none",
                plot.caption = element_text(face = "bold", hjust = 0), text=element_text(size=30)) +
          # lims( y = c(0, max(current_index$alphadiv_values)+0.01) ) +
          geom_signif(
            comparisons = list(c(group1, group2)), 
            test = "wilcox.test",
            map_signif_level = TRUE,
            y_position = max(current_index$alphadiv_values)+0.01)
        
        
        tiff(paste(output,"/alpha/alpha_diversity_",index, "_", group1, "_",group2,".tiff", sep=""), width=600, height=600)
        print(alpha_plot)
        dev.off()

        used_groups = append(used_groups, group1)
        

    }}}
  }
    
}




row.names(data_cols) <- metadata$Taxa

## alpha diversity at the Genus level
# df <- data_cols |> t() |> as.data.frame()
# # df$Species = row.names(df)
# df$Genera <- sapply(str_split(row.names(df), ";"), `[`, 1)
# df <- aggregate(df, . ~ Genera, FUN = sum)
# df <- df[df$Genera %in% df$Genera[grep("g__", df$Genera)], ]
# row.names(df) = df$Genera
# df$Genera = NULL
# df <- df |> t() |> as.data.frame()
# df$Taxa = row.names(df)
# 
# alpha_diversity(df, metadata, "Group_Condition")
# 
# 
# tongue_metadata <- metadata[metadata$Group_Site == "Xerostomic-T" | metadata$Group_Site == "Control-T",]
# 
# df[ c(tongue_metadata$Taxa), ] %>% View()
# 
# alpha_diversity(df[ c(tongue_metadata$Taxa), ], tongue_metadata, "Group_Site")


# alpha_diversity(data_cols, metadata, "Group_Condition")


# tongue_metadata <- metadata[metadata$Group_Site == "Xerostomic-T" | metadata$Group_Site == "Control-T",]

# data_cols[ c(tongue_metadata$Taxa), ] %>% View()

# alpha_diversity(data_cols[ c(tongue_metadata$Taxa), ], tongue_metadata, "Group_Site")

alpha_diversity(data_cols, metadata, "Group")

