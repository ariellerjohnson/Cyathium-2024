#!/bin/bash

#SBATCH --job-name=fastp.sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arj66@cornell.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30gb
#SBATCH --time=999:00:00
#SBATCH --output=fastp.joblog

#file location and project prefix
RNAseq_raw_dir=/home/FrankLab/Arielle/Euphorbia_cyathium
project_name=cyathium
underscore_field=8

#print the date and time
date

#make directory for this step
mkdir fastp_trim
cd fastp_trim

#do a for loop to get the files
RNAseq_files=( $(ls $RNAseq_raw_dir/*fastq.gz))

#run fastp with parameters --cut_window_size 5 --cut_mean_quality 20 --length_required 50
for file in ${RNAseq_files[*]}
do
  eval "prefix=$(echo ${file} | cut -d_ -f $underscore_field)"
  eval "my_prefix=${project_name}_${prefix}"
  /programs/fastp-0.23.2/fastp --cut_right --cut_window_size 5 --cut_mean_quality 20 --length_required 50 -i ${file} -o ${my_prefix}.fq -j ${my_prefix}.json -h ${my_prefix}.html
done

#print the date and time again
date
