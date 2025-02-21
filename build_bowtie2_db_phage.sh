#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=bowtie2_index          ## Name of the job
#SBATCH -A katrine_lab                    ## Account to charge
#SBATCH -p standard                       ## Partition/queue name
#SBATCH --ntasks=1                        ## Number of processes to be launched
#SBATCH --mem=4G                          ## Memory request
#SBATCH --cpus-per-task=4                 ## Number of cores per task
#SBATCH --error=bowtie2_index-%J.err      ## Error log file
#SBATCH --output=bowtie2_index-%J.out     ## Output info file

# Load Bowtie2 module
module load bowtie2/2.4.5

# Set the path to your reference genome FASTA file and output directory
REF_FASTA="/pub/amonsiba/Projects/240820_RNAseq.ch2/ref/250113_SBP2phi2_ref/SBP2phi2_terminase_adjusted/SBP2phi2_reordered.fasta"                         # Your genome assembly file in FASTA format
OUT_PREFIX="/pub/amonsiba/Projects/240820_RNAseq.ch2/ref/250113_SBP2phi2_ref/SBP2phi2_terminase_adjusted/SBP2phi2_bowtie_db"                             # Prefix for the Bowtie2 index files

# Step 1: Index the genome assembly with Bowtie2
echo "Indexing genome assembly with Bowtie2..."
bowtie2-build $REF_FASTA $OUT_PREFIX

# Check if the indexing was successful
if [ $? -eq 0 ]; then
  echo "Bowtie2 indexing completed successfully."
else
  echo "Error: Bowtie2 indexing failed."
  exit 1
fi
