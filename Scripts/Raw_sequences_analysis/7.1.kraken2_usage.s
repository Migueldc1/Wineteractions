#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J kraken2_tot
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/kraken2_tot.%j.out
#SBATCH -e slurm/kraken2_tot.%j.err
#SBATCH --mail-user migueldc@ucm.es
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL

for read1 in reads/*R1_001_72nt.fastq.gz; 
do read2=$(echo $read1| sed s/'_R1_'/'_R2_'/); 
sample_name=$(echo $read1| grep -oP '(?<=/).*(?=_S)'); 
echo $read1; 
echo $read2; 
echo $sample_name;

$HOME/kraken2/kraken2 --use-names --threads 16 --db Kraken2/Wine/vino --report kraken2_out/$sample_name.report.txt --paired $read1 $read2 > kraken2_out/$sample_name.kraken ;

done