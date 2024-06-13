#!/bin/bash -l

#SBATCH --job-name=mapDamage
#SBATCH -o mapDamage-%j.out
#SBATCH --time=144:00:00
#SBATCH -p main
#SBATCH -c 4

enable_lmod
module load container_env mapdamage2

ls PIRE2019-Ssp-*noRS.bam | parallel --no-notice -kj40 "crun mapDamage -i{} -r ../../reference.5.5.fasta -l 150" 
