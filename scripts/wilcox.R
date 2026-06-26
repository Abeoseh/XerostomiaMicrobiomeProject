## library packages 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library())


## set working directory
# getwd()
# setwd("X:/OralMedRsch/Abeoseh Flemister/XerostomiaMicrobiomeProject/")
# getwd()


## set output dir
# output="./output/without_men"
output="./output/without_men_10percent_filtered"
# output="./output/without_men_site"


## read file 
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx")


microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()


microbiome_data$Group[ grep("Control", microbiome_data$Group) ] <- "Non-Xerostomia"
microbiome_data$Condition[microbiome_data$Condition == "None"] <- ""
microbiome_data$Group_Condition <- paste(microbiome_data$Group, microbiome_data$Condition)
microbiome_data$Group_Site <- paste(microbiome_data$Group, microbiome_data$Site)
microbiome_data$Group_Condition_Site <- paste(microbiome_data$Group_Condition, "-", microbiome_data$Site, sep = "")
microbiome_data <- relocate(microbiome_data, Group_Condition, Group_Site, Group_Condition_Site)
# head(microbiome_data)

# groups <- unique(microbiome_data$Group)



calculate_wilcox <- function(microbiome_data, column){   
    groups <- unique(microbiome_data[[column]])
    taxa_df <- data.frame()
    ## mann whitney-u on individual taxa
    used_groups = c()
    for (group1 in groups){
      
      for (group2 in groups){
        if (group1 != group2 & !(group2 %in% used_groups) ){
          
          for (taxa_index in 9:length(microbiome_data)){
            
            taxa_at_index = colnames(microbiome_data)[taxa_index]
            # print(taxa_at_index)
            current_taxa <- microbiome_data[,c(column,taxa_at_index)] ## select the current taxa and the groups column
            
            
          
            taxa1 <- current_taxa[current_taxa[[column]] == group1,2] %>% unlist() ## Control
            taxa2 <- current_taxa[current_taxa[[column]] == group2,2] %>% unlist() ## Xerostomia
            mean_group1 <- mean(taxa1)
            mean_group2 <- mean(taxa2)
            ## taxa1 + 1 because a lot of the values are 0.. positive log2fold change means the bacteria was found in greater quantities in Xerostomia
            # logchange <- log2(abs( (mean(taxa1 + 1) - mean(taxa2 + 1)) / mean(taxa2 + 1) ) ) 
            logchange <- mean(log2(taxa1 + 1)) - mean(log2(taxa2 + 1)) ## mathemathically the same as log2(geo_mean(condition1 + 1)/geo_mean(condition2 + 1))
            
            pvalue <- wilcox.test(taxa1, taxa2, alternative = "two.sided", paired = FALSE)$p.value #
            # pvalue_corrected <- pvalue * (length(microbiome_data)- 5) ##bonferonni correction
            
            taxa_df <- rbind(taxa_df, data.frame(bacteria = colnames(current_taxa)[2], group1 = group1, 
                                                 group2 = group2, group1_mean = mean_group1, group2_mean = mean_group2, log2foldchange = logchange, pvalue = pvalue, pvalue_corrected=NA)  )
            
            
  
          }
          
          taxa_df$pvalue_corrected = p.adjust(taxa_df$pvalue, method = "BH")
          df = taxa_df    
          
          # df$pvalue_corrected <- p.adjust(taxa_df$pvalue, method = "BH")
          # print(head(df))
          # return(filter(taxa_df, pvalue_corrected < 0.049))
          # filter(df,group1 == !!group1  & group2 == !!group2) |> print()
          df = filter(df, pvalue_corrected < 0.049 & group1 == !!group1 & group2 == !!group2)
          # print(df)
          ## select the top 10 and bottom 10 logFC
          wilcox_plot = ggplot(rbind(df |> arrange(desc(log2foldchange)) |> slice(1:10), df |> arrange((log2foldchange)) |> slice(1:10)), aes(log2foldchange, bacteria, color = group1_mean - group2_mean > 0)) + 
            geom_point() +
            scale_color_hue(name = column, labels = c(unique(df$group2), unique(df$group1))) +
            geom_vline(xintercept = 0) +
            labs(title = "Top and Bottom 10 log2FC", x = "log2FC", y = "Bacteria")
          
          dir.create(file.path( paste(output,"/wilcox/plots", sep=""),column ), showWarnings = FALSE)
          png(paste(output,"/wilcox/plots/",column,"/log2FC_", unique(df$group1), "_",unique(df$group2),".png", sep=""), width = 700)
          print(wilcox_plot)
          dev.off()
          rm(df)
          }
        
        }
      used_groups = append(used_groups, group1)
      
      

      
    }
    

    # taxa_df$pvalue_corrected <- p.adjust(taxa_df$pvalue, method = "BH")
    write.csv(taxa_df, paste(output,"/wilcox/wilcox_pvalue_corrected_benjamin hochberg_",column,".csv",sep=""), row.names = F)
    write.csv(filter(taxa_df, pvalue_corrected < 0.049), paste(output,"/wilcox/wilcox_pvalue_corrected_benjamin hochberg_",column, "_sigonly.csv",sep=""), row.names = F)
    
}


calculate_wilcox(microbiome_data, colnames(microbiome_data)[1]) ## Group_Condition
calculate_wilcox(microbiome_data, colnames(microbiome_data)[2]) ## Group_Site
calculate_wilcox(microbiome_data, colnames(microbiome_data)[3]) ## Group_Condition_Site
calculate_wilcox(microbiome_data, colnames(microbiome_data)[4]) ## Group
calculate_wilcox(microbiome_data, colnames(microbiome_data)[5]) ## Condition
calculate_wilcox(microbiome_data, colnames(microbiome_data)[6]) ## Site




# taxa_df %>% View()
# filter(taxa_df, pvalue_corrected < 0.049) %>% View() 
 


# df <- read.csv("./output/lefse/lefseresults.csv")
# common_bacteria = intersect(filter(taxa_df, pvalue_corrected < 0.049 & (log2foldchange < 0 | log2foldchange > 0))$bacteria,
#                                    df$features )
# write.csv(filter(taxa_df, bacteria %in% common_bacteria), "./output/wilcox/lefse_wilcox_intersect.csv", row.names =  F)

# ggplot(filter(taxa_df, bacteria %in% common_bacteria), aes())

