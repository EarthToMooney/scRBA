#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --partition=open 
module load gams
outputFolder=\""./$(basename $(pwd))/\""
cd ../
python A1_run_yield_predictions.py outputFolder=$outputFolder
