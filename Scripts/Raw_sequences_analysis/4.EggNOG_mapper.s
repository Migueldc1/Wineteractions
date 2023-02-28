#!/bin/bash
#SBATCH -p long
#SBATCH -n 8
#SBATCH -J eggnog
#SBATCH --mem 12288
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/eggnog.%j.out
#SBATCH -e slurm/eggnog.%j.err

source /home/${USER}/.bashrc
source activate eggnog

for smpl in Alignment/Index/*; 
do sample_name=$(echo $smpl| grep -oP '(?<=/).*'|  grep -oP '(?<=/).*');  
echo $sample_name;

mkdir Annotation/eggnog/$sample_name

emapper.py -m diamond --itype CDS -i Alignment/Index/$sample_name/$sample_name".fa" -o Annotation/eggnog/$sample_name/$sample_name --decorate_gff yes --dmnd_db Databases/eggnog/saccharomycetaceae.dmnd;

gffread Annotation/eggnog/$sample_name/$sample_name".emapper.decorated.gff" -F -T -o $sample_name".gtf";

done
