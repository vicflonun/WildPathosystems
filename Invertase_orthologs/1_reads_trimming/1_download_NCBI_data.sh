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
#module load $MODULE_NAME 

#########################
# Run commands from here:


# Download Illumina reads for Zt549 
#BioProject PRJNA1162695; BioSample SAMN43819765; Experiment SRX26210164; Run (reads) SRR30809458 

mkdir -p NCBI_data
cd NCBI_data

fasterq-dump SRR30809458 --split-files --progress

gzip *fastq


