#!/bin/bash

#SBATCH --job-name=md5sum       ## Name of the job.
#SBATCH -A katrine_lab      ## CHANGE account to charge
#SBATCH -p standard               ## partition name
#SBATCH --nodes=1             ## (-N) number of nodes to use
#SBATCH --ntasks=1            ## (-n) number of tasks to launch
#SBATCH --array=1-40          ## array
#SBATCH --cpus-per-task=1     ## number of cores the job needs
#SBATCH --error=slurm-%J.err  ## error log file
#SBATCH --output=slurm-%J.out ## output log file

# Set the base directory (adjust if necessary)
FASTQ="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/rawreads"
OUT="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/qc"

#Make a temp file
temp=$(basename -s .fastq.gz $FASTQ/*.fastq.gz | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

md5sum $FASTQ/${prefix}.fastq.gz >> $OUT/"240820_md5sumfile.txt"
