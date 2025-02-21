#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=clean_rnaseqreads_1    ## Name of the job.
#SBATCH -A katrine_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1-20   ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=10    ## number of cores the job needs
#SBATCH --mem=6G
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file

#load modules
module load fastqc/0.11.9
module load bbmap/38.87

# Set directory paths
RAW="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/rawreads"
CLEAN="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/cleanreads"
PREQC="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/qc/preclean"
POSTQC="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/qc/postclean"

#Make a temp file
temp=$(basename -s _R1_001.fastq.gz $RAW/*_R1_001.fastq.gz | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# Step 1: Initial Quality control
echo "Running Initial FastQC on raw reads..."
fastqc $RAW/${prefix}_R1_001.fastq.gz $RAW/${prefix}_R2_001.fastq.gz \
	-o $PREQC/

# Step 2: Trimming and quality filtering
echo "Trimming adapters and quality filtering with BBDuk..."
bbduk.sh in1=$RAW/${prefix}_R1_001.fastq.gz in2=$RAW/${prefix}_R2_001.fastq.gz \
         ref=adapters,phix \
	 ktrim=r \
	 mink=11 \
	 hdist=1 \
         qtrim=rl \
         ftl=15 \
         ftr=149 \
         trimq=30 \
         minlen=10 \
         out=$CLEAN/${prefix}_R1_trim.fastq.gz \
         out2=$CLEAN/${prefix}_R2_trim.fastq.gz \
         stats=$CLEAN/${prefix}_bbduk_stats1.txt \
         refstats=$CLEAN/${prefix}_bbduk_ref_stats1.txt \
         tpe \
         tbo \
         threads=$SLURM_CPUS_PER_TASK

echo "RNAseq Cleaning pipeline #1  completed!"


