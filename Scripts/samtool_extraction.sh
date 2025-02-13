#!/bin/bash
#SBATCH --job-name=samtool_extraction
#SBATCH -A katrine_lab     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --output=logs/samtoolextraction_%A_%a.out
#SBATCH --error=logs/samtoolextraction_%A_%a.err

# Define the directory containing your SAM files
SAM="/pub/amonsiba/Projects/240820_RNAseq.ch2/results/coverage/Phage/txt_files"

output_file="$SAM/coverage_summary_phage.csv"

# Write header to output file
echo "Sample,Number_of_Reads,Percent_Covered,Mean_Coverage" > "$output_file"

# Loop through each coverage file
for file in "$SAM"/*.coverage.txt; do
    # Extract sample name from file name
    sample_name=$(basename "$file" .coverage.txt)

    # Extract metrics
    num_reads=$(grep "Number of reads:" "$file" | awk '{print $NF}')
    percent_covered=$(grep "Percent covered:" "$file" |  awk '{print $6}' | grep "%")
    mean_coverage=$(grep "Mean coverage:" "$file" |  awk '{print $6}' | grep "x")

    # Append the results to the output file
    echo "$sample_name,$num_reads,$percent_covered,$mean_coverage" >> "$output_file"
done

echo "Data extraction complete. Results saved to $output_file."

