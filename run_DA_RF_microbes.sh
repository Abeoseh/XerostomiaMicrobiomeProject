#!/bin/bash

module load R/4.5.1

Rscript ./scripts/plot_DA_RF_microbes.R \
		"csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx" \
		"./output/without_men_saliva_10percent_filtered/roc/auc_df.csv" \
		./output/without_men_saliva_10percent_filtered/ROC_fisher

Rscript ./scripts/plot_DA_RF_microbes.R \
		"csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx" \
		"./output/without_men_10percent_filtered/roc/auc_df.csv" \
		./output/without_men_10percent_filtered/ROC_fisher
