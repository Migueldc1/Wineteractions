#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J asmbl
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/transabyss.%j.out
#SBATCH -e slurm/transabyss.%j.err

module swap openmpi openmpi3/3.1.0

source /home/${USER}/.bashrc
source activate transabyss

for read1 in Reads/non-rRNA/*_R1*; 
do read2=$(echo $read1| sed s/'_R1'/'_R2'/);  
sample_name=$(echo $read1| grep -oP '(?<=/).*' | grep -oP '(?<=/).*(?=_R)'); 
echo $sample_name;  
echo $read1;  
echo $read2;

mkdir Assembly/transabyss/$sample_name

transabyss --pe $read1 $read2 --outdir Assembly/transabyss/$sample_name/kmer21 --kmer 21;
transabyss --pe $read1 $read2 --outdir Assembly/transabyss/$sample_name/kmer29 --kmer 29;
transabyss --pe $read1 $read2 --outdir Assembly/transabyss/$sample_name/kmer39 --kmer 39;
transabyss --pe $read1 $read2 --outdir Assembly/transabyss/$sample_name/kmer59 --kmer 59;

transabyss-merge Assembly/transabyss/$sample_name/kmer21/transabyss-final.fa Assembly/transabyss/$sample_name/kmer29/transabyss-final.fa \
Assembly/transabyss/$sample_name/kmer39/transabyss-final.fa Assembly/transabyss/$sample_name/kmer59/transabyss-final.fa --mink 21 --maxk 59 \
--out Assembly/$sample_name.fa;

done