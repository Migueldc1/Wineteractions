#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J bwa_idx
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/bwa_idx.%j.out
#SBATCH -e slurm/bwa_idx.%j.err

for asmbl in Assembly/*.fa; 
do read1=$(echo $asmbl"_R1.fq.gz");
read2=$(echo $asmbl"_R2.fq.gz");
sample_name=$(echo $asmbl| grep -oP '(?<=/).*' | sed s/'.fa'/''/);  
echo $sample_name;  
echo $read1;  
echo $read2;

mkdir Alignment/Index/$sample_name
mv $asmbl Alignment/Index/$sample_name/$sample_name".fa";

bwa index Alignment/Index/$sample_name/$sample_name".fa";

done