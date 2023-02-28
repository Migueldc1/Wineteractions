#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J assembly-stats
#SBATCH --mem 4096
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/ass-stats.%j.out
#SBATCH -e slurm/ass-stats.%j.err

source /home/${USER}/.bashrc
source activate busco

for asmbl in Assembly/*;
do sample_name=$(echo $asmbl| grep -oP '(?<=/).*' | sed s/'.fa'/''/); 
echo $sample_name;  
echo $asmbl;

assembly-stats -t $asmbl > Assembly/stats/$sample_name".txt"

done