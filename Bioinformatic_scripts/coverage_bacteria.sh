#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=coverage_bacteria_analysis          ## Name of the job
#SBATCH -A katrine_lab                    ## Account to charge
#SBATCH -p standard                           ## Partition/queue name
#SBATCH --ntasks=1                            ## Number of processes to be launched
#SBATCH --array=1-20   ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --mem=6G                              ## Memory request
#SBATCH --cpus-per-task=4                     ## Number of cores per task the job needs
#SBATCH --error=slurm-%J.err                  ## Error log file
#SBATCH --output=slurm-%J.out                 ## Output info file

# Load necessary modules
module load samtools/1.15.1
module load bowtie2/2.4.5

# Set paths to your files
REF="/pub/amonsiba/Projects/240820_RNAseq.ch2/ref/250203_B28B_ref/trycyler_assembly"     # Path to your genome assembly in FASTA format
READS="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/cleanreads" # Path to your FASTQ file (R2), if paired-end
BAM="/pub/amonsiba/Projects/240820_RNAseq.ch2/results/coverage/250203_Coverage_Bacteria/bam"
TXT="/pub/amonsiba/Projects/240820_RNAseq.ch2/results/coverage/250203_Coverage_Bacteria/txt_files"

#Make a temp file
temp=$(basename -s _R1_clean.fastq.gz $READS/*_R1_clean.fastq.gz | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# Step 1: Align reads to the genome using Bowtie2
echo "Aligning reads with Bowtie2..."
bowtie2 -x $REF/B28B_bowtie_db \
 -1 $READS/${prefix}_R1_clean.fastq.gz \
 -2 $READS/${prefix}_R2_clean.fastq.gz \
 --no-unal \
 -p $SLURM_CPUS_PER_TASK |\
samtools view -bS - |\
samtools sort - -o $BAM/${prefix}.sorted.bam
samtools index $BAM/${prefix}.sorted.bam

# Step 2: Calculate coverage statistics using bedtools
echo "Calculating coverage..."
#samtools depth – computes the read depth at each position or region
samtools depth -a $BAM/${prefix}.sorted.bam -o $TXT/${prefix}.depth.txt

#samtools coverage – produces a histogram or table of coverage per chromosome
samtools coverage $BAM/${prefix}.sorted.bam -m -o $TXT/${prefix}.coverage.txt

echo "Coverage analysis completed. Results are in $TXT."

