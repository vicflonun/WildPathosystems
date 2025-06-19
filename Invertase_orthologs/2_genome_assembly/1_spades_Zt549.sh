#!/bin/bash
 
#SBATCH --job-name=spades #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=110G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this fileE
#SBATCH --partition=global #Request a specific partition for the resource allocation.


##########################
# Module load:
module load python/3.7.1

#########################
# Run commands from here:


# Define variables:
workdir="/home/dalsasso/genome_assembly/spades/Zt549"
reads_pe1="/home/dalsasso/data/Zt549_Victor/NCBI_data_cut/SRR30809458_R1.paired.cut.fastq.gz"
reads_pe2="/home/dalsasso/data/Zt549_Victor/NCBI_data_cut/SRR30809458_R2.paired.cut.fastq.gz"


# Run SPAdes:
mkdir -p "$workdir"

/data/biosoftware/spades/SPAdes-3.15.0-Linux/bin/spades.py \
    -t 15 \
    -m 100 \
    --careful \
    -k 21,33,55,67,99,127 \
    -o "$workdir" \
    -1 "$reads_pe1" \
    -2 "$reads_pe2"


# Filter contigs >= 200 bp:
seqtk seq -L 200 "${workdir}/scaffolds.fasta" > "${workdir}/${species}.genome.fasta"



