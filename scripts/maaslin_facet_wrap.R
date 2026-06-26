suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

## read in microbiome data
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()
microbiome_data$Condition[microbiome_data$Condition == "None"] <- "Control"

colnames(microbiome_data)[1:5]
colnames(microbiome_data)[1:6]


combine_maaslin <- function(microbiome_data, metadata_col, maaslin_file_path, output_file_path, group1 = NULL, group2 = NULL){
    if(!is.null(group1)){
      microbiome_data = microbiome_data[microbiome_data[[metadata_col]] == group1 | microbiome_data[[metadata_col]] == group2,]

    }
    
    ## get only the counts
    countdata=microbiome_data[,c(metadata_col, colnames(microbiome_data)[6:length(microbiome_data)]) ] %>% t() %>% as.data.frame()
    row.names(countdata) <- gsub(";|-|\\*", ".", row.names(countdata))
    countdata$feature = row.names(countdata)
    

    
    ## get the Group column and match it to the columns (v1,...,v135) of countdata
    group <- countdata |> t() |> as.data.frame() |> select({{metadata_col}})
    group$bacteria = row.names(group)
    group <- group[1:((dim(group)[1])-1), ]
    

    
    ## read in the maaslin2 df
    df <- read.csv(maaslin_file_path, sep = "\t")
    ## remove group column and merge
    df$Group <- NULL
    df2 <- merge(df, countdata)
    

    
    ## make the df long and add the Group column from microbiome_data
    df3 <- df2 |> 
      pivot_longer(-c(feature, metadata, value, coef, stderr, N, N.not.0, pval, qval), 
                   names_to = "bacteria", values_to = "count") |>
      merge(group, by = "bacteria")

    
    metadata_col_var <- sym(metadata_col)
    ## plot the results
    png(output_file_path, height = 5000, width = 5000)
    plot = ggplot(df3, aes(x = !!metadata_col_var, y = as.numeric(count))) + 
      geom_boxplot(aes(fill = !!metadata_col_var)) + geom_point(aes(color = !!metadata_col_var)) +
      # geom_boxplot(aes(fill = !!metadata_col_var), outliers = FALSE) +
      facet_wrap( ~ feature, ncol = 5, scales = "free") +
      theme(strip.text = element_text(size=25)) +
      # geom_label(aes(label = paste("FDR", qval) )) +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2")
    
    print(plot)
    dev.off()

}



## Xerostomic vs Control
combine_maaslin(microbiome_data, "Group", "./output/maaslin/GroupXerostomicControl/significant_results.tsv", "output/maaslin/GroupXerostomicControl/GroupXerostomicControl.png")

## Medication vs Control
combine_maaslin(microbiome_data, "Condition", 
                "./output/maaslin/ConditionMedicationControl/significant_results.tsv", 
                "output/maaslin/ConditionMedicationControl/ConditionMedicationControl.png", "Medication", "Control")

## Autoimmune vs Medication
combine_maaslin(microbiome_data, "Condition", 
                "./output/maaslin/ConditionAutoimmuneMedication/significant_results.tsv", 
                "output/maaslin/ConditionAutoimmuneMedication/ConditionAutoimmuneMedication.png", "Autoimmune", "Medication")


## Autoimmune vs Control
combine_maaslin(microbiome_data, "Condition", 
                "./output/maaslin/ConditionAutoimmuneControl/significant_results.tsv", 
                "output/maaslin/ConditionAutoimmuneControl/ConditionAutoimmuneControl.png", "Autoimmune", "Control")


