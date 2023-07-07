#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J bracken_tot
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/bracken_tot.%j.out
#SBATCH -e slurm/bracken_tot.%j.err
#SBATCH --mail-user migueldc@ucm.es
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL

for report in kraken2_out/*.report.txt; 
do outp=$(echo $report| grep -oP '(?<=/).*(?=.report)'); 
echo $outp;

$HOME/Bracken-2.5/bracken -d Kraken2/Wine/vino -i $report -r 72 -l S -o bracken_out/$outp.bracken;

done