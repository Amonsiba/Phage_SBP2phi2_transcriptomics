#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=merge_bam
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --array=1-20         ## Adjust array size to match the number of unique prefixes
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --output=bammerge_%A_%a.out
#SBATCH --error=bammerge_%A_%a.err

# Load necessary modules if required
module load bowtie2/2.5.1 
module load samtools/1.15.1

# Define the directory containing your BAM files
BAM="/pub/amonsiba/Projects/240820_RNAseq.ch2/results/coverage/250113_Coverage_Phage/terminase_adjusted_fasta/bam"

# Step 1: Get unique prefixes for samples (e.g., B28B_t0, Ph_t10, etc.)
prefixes=$(ls $BAM/*.sorted.bam | sed -E 's/.*\/(.*_t[0-9]+)_rep[0-9]_S[0-9]+\.sorted\.bam/\1/' | sort -u)

# Convert the list of prefixes into an array
prefix_array=($prefixes)

# Get the prefix for the current task
prefix=${prefix_array[$SLURM_ARRAY_TASK_ID-1]}

# Step 2: Find BAM files for rep1 and rep2 for the current prefix
bam_files=$(ls $BAM/${prefix}_rep*.sorted.bam)

# Step 3: Merge BAM files
merged_bam=$BAM/${prefix}_merged.bam
echo "Merging BAM files for prefix: $prefix"
samtools merge $merged_bam $bam_files

# Step 4: Index the merged BAM file
echo "Indexing merged BAM file: $merged_bam"
samtools index $merged_bam

echo "Finished processing prefix: $prefix"
