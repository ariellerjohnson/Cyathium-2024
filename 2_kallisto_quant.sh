#!/bin/bash

#SBATCH --job-name=kallisto_quant.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30gb
#SBATCH --time=999:00:00
#SBATCH --output=kallisto_quant.joblog

#file locations and read length parameters
RNADIR=/workdir/arj66/Cyathium_RNAseq/fastp_trim #use the location of the trimmed files
REFERENCE=/workdir/arj66/Cyathium_RNAseq/Euphorbia_peplus.codingseq
MEAN_LENGTH=275
SD=86

#print the date and time
date

#make directory for this step
mkdir kallisto_count
cd kallisto_count

#make the index file from the reference transcriptome
/programs/kallisto/kallisto index -i transcripts.idx $REFERENCE

#use a for loop to run quantification on all of the RNAseq files
RNA_FILES=$RNADIR/*.fq

for RNA_FILE in ${RNA_FILES[*]}
do
  eval "OUTFOLDER_NAME=$(echo ${RNA_FILE##*/} | cut -d. -f 1)"
  echo $OUTFOLDER_NAME
  /programs/kallisto/kallisto quant --single -i transcripts.idx -o $OUTFOLDER_NAME -b 100 -t 10 -l $MEAN_LENGTH -s $SD $RNA_FILE
done

#print the date and time again
date
