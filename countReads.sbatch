#!/bin/bash

#SBATCH --job-name=cnt
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1


module load parallel

# Count reads
ls *2019*[FR].fq.gz | parallel "echo -n {}' ' && zgrep -c '^\+$' {}" | sort -g >> Untrimmed_NumReads.txt
