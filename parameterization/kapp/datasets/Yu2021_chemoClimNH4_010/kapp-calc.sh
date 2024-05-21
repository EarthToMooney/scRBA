#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --partition=open
python3 A1_process_data.py
python3 B1_enz_from_proteome.py
python3 B2_min_flux_violation.py
python3 B3_min_flux_sum.py
python3 C1_calculate_kapp.py
