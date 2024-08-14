#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --partition=open 
module load gams
gams fva-prosyn-maxonly.gms lo=0 --mu=0.38
