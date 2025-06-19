#!/bin/bash

#BATCH --job-name=JOBNAME #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=10G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


##########################
# Module load:
module load  python/3.6.7

#########################
# Run commands from here:


# Define the paths
reads_dir="./NCBI_data"        
output_dir="./NCBI_data_cut" 

# Adapter sequences
adapter_fwd="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter_rev="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Parameters for cutadapt
min_length=50
max_n=2
quality_threshold=30
threads=6


mkdir -p "$output_dir"

for forward_read in "$reads_dir"/*_1.fastq.gz; do
    reverse_read=${forward_read/_1.fastq.gz/_2.fastq.gz}
    base_name=$(basename "$forward_read" _1.fastq.gz)

    # Run cutadapt
    cutadapt \
        -a "$adapter_fwd" -b "$adapter_rev" \
        -m $min_length --max-n $max_n -q $quality_threshold -j $threads \
        -o "$output_dir/${base_name}_R1.paired.cut.fastq.gz" \
        -p "$output_dir/${base_name}_R2.paired.cut.fastq.gz" \
        "$forward_read" "$reverse_read" \
        > "$output_dir/${base_name}_cutadapt_summary.txt"
done

