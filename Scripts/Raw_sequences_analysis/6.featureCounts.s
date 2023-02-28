#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J count
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/count.%j.out
#SBATCH -e slurm/count.%j.err

for smpl in Annotation/eggnog/*; 
do sample_name=$(echo $smpl| grep -oP '(?<=/).*'|  grep -oP '(?<=/).*');  
echo $sample_name;

featureCounts -a Annotation/eggnog/$sample_name/$sample_name".gtf" -o Count/$sample_name".txt" Alignment/Output/$sample_name/$sample_name".bam" -t transcript -g em_target;

done