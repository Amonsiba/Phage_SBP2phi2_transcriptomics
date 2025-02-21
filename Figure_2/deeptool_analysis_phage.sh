#!/bin/bash
#SBATCH --job-name=bam_to_bigwig
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --array=1-10
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --output=bam_to_bigwig_%A_%a.out
#SBATCH --error=bam_to_bigwig_%A_%a.err

# Load necessary modules if required
module load deeptools/3.5.1
#source activate deeptools

# Define the directory containing your SAM files
BAM="/pub/amonsiba/Projects/240820_RNAseq.ch2/results/coverage/250113_Coverage_Phage/terminase_adjusted_fasta/bam"
BIGWIG="/pub/amonsiba/Projects/240820_RNAseq.ch2/results/coverage/250113_Coverage_Phage/terminase_adjusted_fasta/bigwig"

#Make a temp file
temp=$(basename -s _merged.bam $BAM/*_merged.bam | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

#Generate BigWig file for coverage analysis
echo "Generating BigWig file $BIGWIG_FILE for coverage analysis..."
bamCoverage --bam $BAM/${prefix}_merged.bam -o $BIGWIG/${prefix}.SeqDepthNorm.bw \
   --binSize 1 \
   --normalizeUsing BPM

echo "Completed processing $SAM_FILE"

