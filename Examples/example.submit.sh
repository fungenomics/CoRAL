#!/bin/sh
#SBATCH --job-name=CoRAL
#SBATCH --account=account 
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=60GB 

# apptainer image
image=path/to/aptainer.sif

# snakefile 
snakefile=path/to/snakefile.master

# config 
config=<path to config file>

# unlock directory in case of previous errors
apptainer exec --contain --cleanenv --pwd "$PWD" $image snakemake -s ${snakefile} --configfile ${snakefile} --unlock 

# run CoRAL  
apptainer exec --contain --cleanenv --pwd "$PWD" $image snakemake -s ${snakefile} --configfile ${config}  --cores 5
