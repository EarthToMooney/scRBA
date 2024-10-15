#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --mem=1GB 
#SBATCH --time=10:00:00 
#SBATCH --partition=open
python copy_from_template.py && python A1_process_data.py && python B1_enz_from_proteome.py && python B2_min_flux_violation.py && python C1_calculate_kapp.py
