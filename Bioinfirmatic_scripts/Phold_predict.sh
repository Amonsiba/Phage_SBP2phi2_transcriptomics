#!/bin/bash

#SBATCH --job-name=phold_predict          ## Job name
#SBATCH -A katrine_lab                          ## Account to charge
#SBATCH -p standard                             ## Partition/queue name
#SBATCH --array=1-8			   	## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --nodes=1                               ## Number of nodes
#SBATCH --ntasks=1                              ## Number of tasks
#SBATCH --cpus-per-task=8                       ## Number of cores (adjust as needed)
#SBATCH --error=slurm-%J.err                    ## Error log file
#SBATCH --output=slurm-%J.out                   ## Output log file

# Define paths and parameters
PHOLD_INPUT="/pub/amonsiba/Projects/250108_PhageRef/Phold_files/genbank_files"   	# Path to the phage FASTA file
PHOLD_OUTPUT="/pub/amonsiba/Projects/250108_PhageRef/Phold_files/Output"      # Output directory for Pharokka results
DB_DIR="/dfs7/whitesonlab/amonsiba/databases/Phold_database"                           # Path to Pharokka database directory
THREADS=8                                                               # Number of threads to use

#Make a temp file
temp=$(basename -s .gbk $PHOLD_INPUT/*.gbk | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`


#Module load phold
echo "Activating pholdENV environment"
eval "$(conda shell.bash hook)"
conda activate /data/homezvol0/amonsiba/.conda/envs/pholdENV

# Run Phold Predict - Step 1 
phold predict -i $PHOLD_INPUT/${prefix}.gbk -o $PHOLD_OUTPUT/Predict/${prefix} -d $DB_DIR -t 8

# Explanation of options:
# -i : Path to input phage GenBank file
# -o : Output directory for Phold results (step1) 
# -t : Number of threads

echo "Phold predict completed for the phage genome."

