#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J eggnog_db
#SBATCH --mem 65536
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/eggnog_db.%j.out
#SBATCH -e slurm/eggnog_db.%j.err

source /home/${USER}/.bashrc
source activate eggnog

mkdir Databases/eggnog
download_eggnog_data.py --data_dir Databases/eggnog -y
