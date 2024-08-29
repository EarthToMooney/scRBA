# Model-specific settings for kapp calculations, as well as steps that must be run each time
import pandas as pd
import numpy as np
import requests
import json
import os,shutil
from gsm_custom_functions import *
from simulate import RBA_result

pycore_path = '../../../../pycore/'
import sys
sys.path.append(pycore_path)
from utils import metabolites_dict_from_reaction_equation_RBA
from utils import extract_details_from_rxnid
vmax = 1e3 # max flux in either direction

# Load path
path_gen = '../../../../build_model/'
path_gams = '../../../../GAMS/'

prot_path = path_gen + 'input/PROTEIN_stoich_curation.xlsx'
# prot_path = path_gen + 'input/PROTEIN_stoich_curation.csv'
unknown_prot_path = path_gen + 'input/PROTEIN_stoich_curation_unknown.xlsx'
model_xlsx_path = path_gams + 'model/RBA_stoichiometry.xlsx'
ribonuc_path = path_gen + 'input/RIBOSOME_nucleus.xlsx'
ribomito_path = path_gen + 'input/RIBOSOME_mitochondria.xlsx'
gsm_rxn_ids_path = path_gams + 'model/GSM_rxn_ids.txt'
sij_path = path_gams + 'model/RBA_sij.txt'
aa_mapping_path = path_gen + 'input/PROTEIN_amino_acid_map.txt'
aa_map = pd.read_csv(aa_mapping_path, sep='\t')
aa_dict = dict(zip(aa_map['aa_abbv'], aa_map['MW']))

flux_path = '../raw_data_files/SC_MFA_flux_CENPK_batch_Rabinowitz.xlsx'
col_LB = 'mfaLB'
col_UB = 'mfaUB'

# must match the growth rate in other files
mu = 0.391-0.024
# protein fraction (disable by uncommenting "ptot = 1" unless composition varies w/ growth rate)
ptot = 0.5551

# Use only name and abundance cols
df_raw = pd.read_excel('../raw_data_files/Shen2024_41589_2024_1571_MOESM3_ESM.xlsx',
                         sheet_name='Table 10a. abs_prot_SC_CENPK')
# df_raw = pd.read_excel('../raw_data_files/Rekena_Datasets.xlsx',
#                          sheet_name='S2 Dataset Final', usecols=[0,1,2,3,4,5,6,17,18,19,20,21,22])
df_raw.index = df_raw['geneID'].to_list() # name of protein
#cols_data = ['replicate 1','replicate 2','replicate 3','replicate 4'] # where protein abundance data is stored
# If encountering issues with averaging multiple replicates, could try using only one column in case the average results in an infeasible AA composition
cols_data = ['replicate 1 normalized to robustly detected frac (93.6%)'] # where protein abundance data is stored
uniprot_col = 'Entry' # set to '' if no column with accession names provided
data_uses_biomass_mass_fraction = False # True if using units like g protein/gDW, False if using g/g protein

dummy_protein = {'id':'PROSYN-PROTDUMMY','AA abundances':dict()}

path_data = './Shen2024_batch_glc.xlsx'

biom_id = 'BIOSYN-BIODILAERO-NOGAM'