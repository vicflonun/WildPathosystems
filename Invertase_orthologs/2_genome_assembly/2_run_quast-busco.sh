#!/bin/bash

#BATCH --job-name=quast #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=20G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=dalsasso@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.




##########################
# Load modules:
module load java/x64/8u121
module load perl/5.24.1
module load python/3.6.0

##########################
# Set variables:

species="Zt549"
base_dir="/home/dalsasso/genome_assembly/spades/${species}"
assembly_file="${base_dir}/${species}.genome.fasta"
quast_output="${base_dir}/quast_${species}"
busco_output="${base_dir}/busco_${species}"
threads="10"

reads_pe1="/home/dalsasso/data/Zt549_Victor/NCBI_data_cut/SRR30809458_R1.paired.cut.fastq.gz"
reads_pe2="/home/dalsasso/data/Zt549_Victor/NCBI_data_cut/SRR30809458_R2.paired.cut.fastq.gz"


# Run QUAST:
quast.py -o "$quast_output" --fungus --fragmented "$assembly_file" \
    --pe1 "$reads_pe1" --pe2 "$reads_pe2" \
    --circos -t "$threads"


# Run BUSCO:
#source ~/.bashrc
conda activate busco

busco -i "$assembly_file" \
    -f -o "$busco_output" \
    -l dothideomycetes_odb10 \
    -m genome -c "$threads"

conda deactivate



