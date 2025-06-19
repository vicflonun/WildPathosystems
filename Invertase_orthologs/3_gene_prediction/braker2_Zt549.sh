#!/bin/bash
 
#SBATCH --job-name=BRAKER2 #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=32G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --partition=global #Request a specific partition for the resource allocation.


################################
# Module load and env preparation

# Get gm_key for GeneMark at: http://topaz.gatech.edu/GeneMark/license_download.cgi
# Each key is valid for limited period and should be replaced after it

#cp /data/biosoftware/braker2/deps/gm_key_64 $HOME/.gm_key
module load perl/5.26.1
source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh
conda activate
module load java/x64/8u121
export PATH=/data/biosoftware/braker2/BRAKER-2.1.6/scripts:$PATH

AUGUSTUS_BASE=/data/biosoftware/braker2/deps/Augustus

export AUGUSTUS_CONFIG_PATH=$AUGUSTUS_BASE/config/
export AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_BASE/scripts
export AUGUSTUS_BIN_PATH=$AUGUSTUS_BASE/bin
export GENEMARK_PATH=/data/biosoftware/braker2/deps/gmes_linux_64/
export BAMTOOLS_PATH=/data/biosoftware/braker2/deps/bamtools/bin/usr/local/bin/
export DIAMOND_PATH=/data/biosoftware/braker2/deps/diamond/
export BLAST_PATH=/data/biosoftware/braker2/deps/ncbi-blast-2.11.0+/bin/
export PROTHINT_PATH=/data/biosoftware/braker2/deps/ProtHint/bin/
export SAMTOOLS_PATH=/data/biosoftware/braker2/deps/samtools/
export CDBTOOLS_PATH=/data/biosoftware/braker2/deps/cdbfasta/

################################

# Inputs
species="Zt549"
workdir="/home/dalsasso/genome_assembly/gene_prediction/${species}"
assembly="/home/dalsasso/genome_assembly/spades/Zt549/Zt549.fasta"

# Input protein databases
speciesDB="/home/dalsasso/data/References/Sources/IPO323_JGI_Lapalu_2023/Zymtr1_INRAE_GeneModels_proteins_20230422.aa.fasta"
# Get orthoDB proteins at: http://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/
orthoDB="/home/dalsasso/data/orthoDB_v12/Fungi.fa"

ncores=10


# Prepare protein DBs
mkdir -p "$workdir"
cd "$workdir" || exit

cat "$speciesDB" "$orthoDB" > "${species}_targetDB.fasta"
PROTEIN_database="${species}_targetDB.fasta"

# Run BRAKER2
braker.pl --fungus \
  --cores "$ncores" \
  --species="$species" \
  --genome="$assembly" \
  --prot_seq="$PROTEIN_database" \
  --UTR=off \
  --gff3 \
  --alternatives-from-evidence=false

rm "$PROTEIN_database"
