#!/bin/bash 
#SBATCH --partition=Orion
#SBATCH --job-name=mafft
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

module unload R
module load R/4.3.3

# RandomForest_ROCabove0.7
# RandomForest

#sbatch sbatch_random_forest.sbatch "csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx" \
#	./output/without_men/RandomForest_ROC_MB

#sbatch sbatch_random_forest.sbatch "csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx" \
#	./output/without_men/RandomForest_ROC_MB

# sbatch sbatch_random_forest.sbatch "csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7.xlsx" \
#	./output/without_men/RandomForest_ROCabove0.7

# sbatch sbatch_random_forest.sbatch "csv_files/Xero-Microbiome RA data patient matched AH 07 23 25 females ROC above 0.7 saliva.xlsx" \
#	./output/without_men_saliva/RandomForest_ROCabove0.7


Rscript ./scripts/random_forest.R &> $2/RF.log "$1" "$2"
# Rscript ./scripts/random_forest_MB.R &> $2/RF.log "$1" "$2"

echo ""
echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
