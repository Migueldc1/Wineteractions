#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J busco
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/busco.%j.out
#SBATCH -e slurm/busco.%j.err

source /home/${USER}/.bashrc
source activate busco

for asmbl in Assembly/*;
do sample_name=$(echo $asmbl| grep -oP '(?<=/).*' | sed s/'.fa'/''/); 
echo $sample_name;  
echo $asmbl;

mkdir Assembly/busco/$sample_name

busco -i $asmbl -l fungi_odb10 -o Assembly/busco/$sample_name -m transcriptome -f

done