#!/bin/bash 
#SBATCH --partition=Orion
#SBATCH --job-name=GenusROC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=100000
echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Num Tasks : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo "" cd $SLURM_SUBMIT_DIR

module load R/4.3.3
Rscript ./scripts/roc_genus.R \
		"csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females saliva 10 percent.xlsx" \
		"./output/without_men_saliva_10percent_filtered"


Rscript ./scripts/roc_genus.R \
		"csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females 10 percent.xlsx" \
		"./output/without_men_10percent_filtered"

echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
