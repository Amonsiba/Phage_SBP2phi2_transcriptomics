#!/bin/bash

#SBATCH --job-name=pharokka_annotation          ## Job name
#SBATCH -A katrine_lab                          ## Account to charge
#SBATCH -p standard                             ## Partition/queue name
#SBATCH --array=1-5			   	## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --nodes=1                               ## Number of nodes
#SBATCH --ntasks=1                              ## Number of tasks
#SBATCH --cpus-per-task=8                       ## Number of cores (adjust as needed)
#SBATCH --error=slurm-%J.err                    ## Error log file
#SBATCH --output=slurm-%J.out                   ## Output log file

# Define paths and parameters
PHAGE_FASTA="/pub/amonsiba/Projects/250108_PhageRef/fasta_files/reordered_fasta_files"   	# Path to the phage FASTA file
OUTPUT_DIR="/pub/amonsiba/Projects/250108_PhageRef/Pharokka_files"      # Output directory for Pharokka results
DB_DIR="/dfs7/whitesonlab/amonsiba/databases"                           # Path to Pharokka database directory
THREADS=8                                                               # Number of threads to use

#Make a temp file
temp=$(basename -s .fasta $PHAGE_FASTA/*.fasta | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`


#Module load pharokka
echo "Activating pharokka environment"
eval "$(conda shell.bash hook)"
conda activate /data/homezvol0/amonsiba/.conda/envs/pharokka

# Run Pharokka
pharokka.py -i $PHAGE_FASTA/${prefix}.fasta -o $OUTPUT_DIR/${prefix}/ -d $DB_DIR -t $THREADS 

# Explanation of options:
# -i : Path to input phage FASTA file
# -o : Output directory for Pharokka results
# -d : Path to the Pharokka database directory (previously installed)
# -t : Number of threads

echo "Pharokka annotation completed for the phage genome."

