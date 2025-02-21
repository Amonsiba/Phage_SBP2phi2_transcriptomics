#!/bin/bash
#--------------------------SBATCH settings------

#SBATCH --job-name=clean_rnaseqreads_2    ## Name of the job.
#SBATCH -A katrine_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=16    ## number of cores the job needs
#SBATCH --mem-per-cpu=5G        ## requested memory (6G = max)
#SBATCH --array=1-20   ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --error=slurm-%J.err ## error log file
#SBATCH --output=slurm-%J.out ##output info file

#load modules
module load fastqc/0.11.9
module load bbmap/38.87
module load bowtie2/2.5.1

# Set directory paths
RAW="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/rawreads"
CLEAN="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/cleanreads"
PREQC="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/qc/preclean"
POSTQC="/pub/amonsiba/Projects/240820_RNAseq.ch2/reads/qc/postclean"

#Make a temp file
temp=$(basename -s _R1_001.fastq.gz $RAW/*_R1_001.fastq.gz | sort -u)
#Make prefix file
prefix=`echo "$temp" | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# Step 3: Deduplication
echo "Removing duplicates with dedupe.sh..."
dedupe.sh \
        in=$CLEAN/${prefix}_R1_trim.fastq.gz \
        in2=$CLEAN/${prefix}_R2_trim.fastq.gz \
        out=$CLEAN/${prefix}.dedupe.fastq.gz \
        threads=$SLURM_CPUS_PER_TASK

# Step 4: Split interleaved output
#echo "Splitting interleaved deduplicated sequences..."
reformat.sh in=$CLEAN/${prefix}.dedupe.fastq.gz \
    out1=$CLEAN/${prefix}_R1_clean.fastq.gz out2=$CLEAN/${prefix}_R2_clean.fastq.gz

# Step 5: Final Quality Check
#echo "Running Final FastQC on raw reads..."
fastqc $CLEAN/${prefix}_R1_clean.fastq.gz $CLEAN/${prefix}_R2_clean.fastq.gz \
	-o $POSTQC/

echo "RNAseq Cleaning pipeline #2 completed!"


