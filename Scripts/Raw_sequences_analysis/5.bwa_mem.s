#!/bin/bash
#SBATCH -p long
#SBATCH -n 16
#SBATCH -J bwa_LM
#SBATCH --mem 32768
#SBATCH -t 7-02:30:15   
#SBATCH -o slurm/bwa_LM.%j.out
#SBATCH -e slurm/bwa_LM.%j.err

mkdir Alignment/Output/LM-CONV-18C

bwa mem Alignment/Index/LM-CONV-18C/LM-CONV-18C.fa Reads/non-rRNA/LM-CONV-18C_R1.fq.gz Reads/non-rRNA/LM-CONV-18C_R2.fq.gz > Alignment/Output/LM-CONV-18C/LM-CONV-18C.sam

mkdir Alignment/Output/LM-CONV-Control

bwa mem Alignment/Index/LM-CONV-Control/LM-CONV-Control.fa Reads/non-rRNA/LM-CONV-Control_R1.fq.gz Reads/non-rRNA/LM-CONV-Control_R2.fq.gz > Alignment/Output/LM-CONV-Control/LM-CONV-Control.sam

mkdir Alignment/Output/LM-CONV-NH4

bwa mem Alignment/Index/LM-CONV-NH4/LM-CONV-NH4.fa Reads/non-rRNA/LM-CONV-NH4_R1.fq.gz Reads/non-rRNA/LM-CONV-NH4_R2.fq.gz > Alignment/Output/LM-CONV-NH4/LM-CONV-NH4.sam

mkdir Alignment/Output/LM-CONV-SO2

bwa mem Alignment/Index/LM-CONV-SO2/LM-CONV-SO2.fa Reads/non-rRNA/LM-CONV-SO2_R1.fq.gz Reads/non-rRNA/LM-CONV-SO2_R2.fq.gz > Alignment/Output/LM-CONV-SO2/LM-CONV-SO2.sam

mkdir Alignment/Output/LM-ECO-18C

bwa mem Alignment/Index/LM-ECO-18C/LM-ECO-18C.fa Reads/non-rRNA/LM-ECO-18C_R1.fq.gz Reads/non-rRNA/LM-ECO-18C_R2.fq.gz > Alignment/Output/LM-ECO-18C/LM-ECO-18C.sam

mkdir Alignment/Output/LM-ECO-Control

bwa mem Alignment/Index/LM-ECO-Control/LM-ECO-Control.fa Reads/non-rRNA/LM-ECO-Control_R1.fq.gz Reads/non-rRNA/LM-ECO-Control_R2.fq.gz > Alignment/Output/LM-ECO-Control/LM-ECO-Control.sam

mkdir Alignment/Output/LM-ECO-NH4

bwa mem Alignment/Index/LM-ECO-NH4/LM-ECO-NH4.fa Reads/non-rRNA/LM-ECO-NH4_R1.fq.gz Reads/non-rRNA/LM-ECO-NH4_R2.fq.gz > Alignment/Output/LM-ECO-NH4/LM-ECO-NH4.sam
