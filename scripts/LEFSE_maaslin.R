.libPaths( c( .libPaths(), "~/my_R_libs4.5.1") )
suppressPackageStartupMessages(library(lefser))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
# suppressPackageStartupMessages(library(Maaslin2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library('SummarizedExperiment'))
# suppressPackageStartupMessages(library(openxlsx))



## set working directory
# getwd()
# setwd("C:/Users/brean/Downloads/masters/Atrium/Xerostomia/")
# getwd()

## set output dir
# output="./output/without_men_saliva"
# output="./output/without_men_10percent_filtered"
output="./output/without_men_saliva_10percent_filtered"



## read file 
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females.xlsx")
# microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx")
microbiome_data <- read_excel("csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx")


microbiome_data$Group[ grep("Control", microbiome_data$Group) ] <- "NX"
microbiome_data$Group[ grep("Xerostomic", microbiome_data$Group) ] <- "XS"
microbiome_data <- relocate(microbiome_data, Group, Condition, Site)
microbiome_data <- microbiome_data[,1:length(microbiome_data)-1] %>% as.data.frame()
microbiome_data$Condition[microbiome_data$Condition == "None"] <- "Control"
# microbiome_data[["XOB ID"]] <- sapply(str_split(microbiome_data$Taxa, "\\."), `[`, 1)


# metadata <- read_excel("csv_files/XOB Clinical Data Microbiome + Metatranscriptomics AH 07 23 2025.xlsx") |>
#   select("Clinical Oral Dryness (COD) Score (1-10)", "XOB ID")
# metadata[["XOB ID"]] <- paste(sapply(str_split(metadata[["XOB ID"]], "_."), `[`, 1),sapply(str_split(metadata[["XOB ID"]], "_."), `[`, 2),sep="0")
# microbiome_data <- merge(microbiome_data, metadata, by = "XOB ID") |> relocate(Group, Condition, Site, "Clinical Oral Dryness (COD) Score (1-10)", "XOB ID")
# microbiome_data$`Clinical Oral Dryness (COD) Score (1-10)` = as.factor(microbiome_data$`Clinical Oral Dryness (COD) Score (1-10)`)



## do LEfSE on the all the metadata columns
metadata_cols <- colnames(microbiome_data)[c(1)]
for (column in metadata_cols){
  

  used_groups = c()
  groups <- unique( microbiome_data[[column]] )
  for (group1 in groups){
    
    for (group2 in groups){
      if (group1 != group2 & !(group2 %in% used_groups) ){
        filtered_microbiome_data <- microbiome_data[microbiome_data[[column]] == group1 | microbiome_data[[column]] == group2,]
        
        countdata=filtered_microbiome_data[,c(8:length(filtered_microbiome_data))] %>% t() %>% as.data.frame()
        colnames(countdata) <- filtered_microbiome_data$Taxa
        ## remove bacteria with NA in them 
        countdata = countdata[!(rownames(countdata) %in% rownames(countdata)[grep("__NA", rownames(countdata))]), ]
        ## select the terminal nodes
        terminal_nodes <- get_terminal_nodes(rownames(countdata))
        countdata <- countdata[terminal_nodes,]
        
        # countdata$Sample <- row.names(countdata)
        # countdata <- relocate(countdata, Sample)
        
        #### format metadata ####
        coldata <- filtered_microbiome_data[,1:6]
        colnames(coldata)[4] = "COD"
        row.names(coldata) <- filtered_microbiome_data$Taxa
        
      
          
          
  
        #### Put into format lefser can understand
        radata=as.matrix(countdata)
        # radata=radata[,-c(1)]
        se=SummarizedExperiment(assays=list(counts=radata),colData=as.matrix(coldata))
        assay(se) = (assay(se))*1e6
        # colData(se) |> View()
        
        se
        
        #### Run lefser ####
        lefser_resultsTP=lefser(se,classCol=column,trim.names=FALSE,checkAbundances=FALSE, lda.threshold=2,kruskal.threshold=0.05,wilcox.threshold=0.05)
        # View(lefser_resultsTP)
        
        ## save results
        # write.csv(lefser_resultsTP, "./output/lefse/lefseresults.csv", row.names = F)
        write.csv(lefser_resultsTP, paste(output, "/lefse/csvs/lefseresults",column,group1,group2,"_bold.csv",sep=""), row.names = F)
        #png(paste(output, "/lefse/plots/lefseresults",column,group1,group2,".png",sep=""), height = 900)
        
        lefser_resultsTP$features <- gsub("g__|f__|o__|c__|p__|k__|NA ", "", gsub(";s__|;", " ", lefser_resultsTP$features))
        
        p <- lefserPlot(lefser_resultsTP,trim.names=FALSE)
        

        
        p$layers[[2]]$aes_params$fontface <- "bold"
        p$layers[[2]]$aes_params$size <- 6
        
        tiff(paste(output, "/lefse/plots/lefseresults",column,group1,group2,"_bold.tiff",sep=""), units="in", width=10, height=15, res=300)
	      print(p)
        dev.off()
      
        
        # #### Run MAASLIN2 ####
        # references = c( "Site,S", "COD,0")
        # for (reference in references){
        #   if (reference == "Site,S"){reference_col = "Site,S"}else{reference_col == "COD,0"}
        #   # # reference = "Site,S"
        #   # # reference_col = "Site"
        #   # 
        #   # reference = "COD,0"
        #   # reference_col = "COD"
        #   
        #   if (column == "Condition"){
        #       reference2 = "Condition,Autoimmune"
        #       countdata2=microbiome_data[,c(8:length(microbiome_data))] %>% t() %>% as.data.frame()
        #       colnames(countdata2) <- microbiome_data$Taxa
        #       ## remove bacteria with NA in them
        #       countdata2 = countdata2[!(rownames(countdata2) %in% rownames(countdata2)[grep("__NA", rownames(countdata2))]), ]
        #       ## select the terminal nodes
        #       terminal_nodes2 <- get_terminal_nodes(rownames(countdata2))
        #       countdata2 <- countdata2[terminal_nodes2,]
        # 
        #       #### format metadata ####
        #       coldata2 <- microbiome_data[,1:6]
        #       colnames(coldata2)[4] = "COD"
        #       row.names(coldata2) <- microbiome_data$Taxa
        # 
        #       fit_data2 = Maaslin2(
        #       input_data = countdata2,
        #       input_metadata = coldata2,
        #       output = paste(output, "/maaslin/",column,"_",reference_col,sep=""),
        #       reference = paste(reference,reference2,sep=";"),
        #       fixed_effects = c(column, reference_col))
        #   }else{
        #     fit_data2 = Maaslin2(
        #       input_data = countdata,
        #       input_metadata = coldata,
        #       output = paste(output, "/maaslin/",column,"_",reference_col,sep=""),
        #       reference = reference,
        #       fixed_effects = c(column, reference_col))
        #     }}
        # 
        # 
        used_groups = append(used_groups, group1)

      }}
    print(paste("done with", column,"column"))}
}



