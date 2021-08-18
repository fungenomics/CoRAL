#!/usr/bin/bash
#PBS -N Zaza_Pipeline
#PBS -o logs/err.txt
#PBS -e logs/out.txt
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=3
#PBS -l mem=125G
#PBS -l vmem=125G


cd /project/kleinman/hussein.lakkis/from_hydra/2021_07_15-Pipeline
mkdir -p logs

#conda init bash
module load miniconda/3.8
source ~/.conda_init.sh
module load snakemake/5.32.0
module load hydra/R/4.0.5
module load python/3.6.5_new

snakemake --use-conda --configfile config.yml --cores 3
