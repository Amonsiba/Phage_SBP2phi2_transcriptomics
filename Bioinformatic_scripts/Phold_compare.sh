#!/bin/bash

#SBATCH --job-name=phold_compare 	         ## Job name
#SBATCH -A katrine_lab                          ## Account to charge
#SBATCH -p standard                             ## Partition/queue name
#SBATCH --array=1-8			   	## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --nodes=1                               ## Number of nodes
#SBATCH --ntasks=1                              ## Number of tasks
#SBATCH --cpus-per-task=8                       ## Number of cores (adjust as needed)
#SBATCH --error=slurm-%J.err                    ## Error log file
#SBATCH --output=slurm-%J.out                   ## Output log file

# Define paths and parameters
GBK="/pub/amonsiba/Projects/250108_PhageRef/Phold_files/genbank_files"   	# Path to the phage genbank file
PREDICT_DIR="/pub/amonsiba/Projects/250108_PhageRef/Phold_files/Output/Predict"
PHOLD_OUTPUT="/pub/amonsiba/Projects/250108_PhageRef/Phold_files/Output/Compare"      # Output directory for Pharokka results
DB_DIR="/dfs7/whitesonlab/amonsiba/databases/Phold_database"                            # Path to Phold database directory
THREADS=8                                                               # Number of threads to use

#Make a temp file
temp=$(basename -s .gbk $GBK/*.gbk | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`


#Module load pharokka
echo "Activating pholdENV environment"
eval "$(conda shell.bash hook)"
conda activate /data/homezvol0/amonsiba/.conda/envs/pholdENV

# Run Phold Compare - Step 2
phold compare -i $GBK/${prefix}.gbk --predictions_dir $PREDICT_DIR/${prefix} -o $PHOLD_OUTPUT/${prefix} -d $DB_DIR -t 8 

# Explanation of options:
# -i : Path to input phage GenBank file
# -o : Output directory for Phold results (step2) 
# -t : Number of threads

echo "Phold compare completed for the phage genome."

