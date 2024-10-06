# Model-specific settings for kapp calculations, as well as steps that must be run each time
import pandas as pd
import numpy as np
import requests,json,os,shutil,sys

model_root_path = '../../../../'
pycore_path = model_root_path+'pycore/'
sys.path.append(pycore_path)
# retrieve name of the directory containing the options file
dir_name = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

from utils import metabolites_dict_from_reaction_equation_RBA, extract_details_from_rxnid
from simulate import * 
from gsm_custom_functions import *
vmax = 1e3 # max flux in either direction, before applying nscale
nscale=1e5 # will be used as scale when running GAMS
nscale2=1e5 # 2nd scale, if different one needed; otherwise, set to same as nscale

# Load path
path_gen = model_root_path+'build_model/'
path_gams = model_root_path+'GAMS/'

prot_path = path_gen + 'input/PROTEIN_stoich_curation.xlsx'
unknown_prot_path = path_gen + 'input/PROTEIN_stoich_curation_unknown.xlsx'
model_xlsx_path = path_gams + 'model/RBA_stoichiometry.xlsx'
ribonuc_path = path_gen + 'input/RIBOSOME_nucleus.xlsx'
ribomito_path = path_gen + 'input/RIBOSOME_mitochondria.xlsx'
gsm_rxn_ids_path = path_gams + 'model/GSM_rxn_ids.txt'
sij_path = path_gams + 'model/RBA_sij.txt'
aa_mapping_path = path_gen + 'input/PROTEIN_amino_acid_map.txt'
aa_map = pd.read_csv(aa_mapping_path, sep='\t')
aa_dict = dict(zip(aa_map['aa_abbv'], aa_map['MW']))

# for flux data; keep path blank if none provided
flux_path = ''
col_LB = 'mfaLB'
col_UB = 'mfaUB'

# must match the growth rate in other files
mu = 0.060 - 0.001
# protein fraction (disable by uncommenting "ptot = 1" unless composition varies w/ growth rate)
ptot = 0.201

# protein data
df_raw = pd.read_excel('../raw_data_files/PinheiroEtAl2020_batchXyl_raw.xlsx',sheet_name='Table S7',header=4, usecols=list(range(31,41)))
# make 1st row the header
df_raw.index = df_raw['IDs in model, or uniprot IDs if not found'].to_list() # name of protein
# df_raw.columns = df_raw.iloc[0]
cols_data = ['REF_P1_g_per_g_prot'] # where protein abundance data is stored
uniprot_col = 'all uniprot matches' # set to '' if no column with accession names provided
data_uses_biomass_mass_fraction = False # True if using units like g protein/gDW, False if using g/g protein
use_ribo_data = False # False to ignore ribosome subunit abundance data

# if True, proteins with measured translation rates of 0 can be made at a flux equal to prosynSlackUB, thus discouraging their use where possible. If you want such proteins to be made, you thus may have to specify them in the phenotype.txt file as having fixed values for prosynSlackUB. 
# if False, proteins with measured translation rates of 0 can't be made. This may cause issues if no isozymes exist for such proteins but their rxns are essential, and such proteins may still be present but unmeasured due to measurement errors.
allow_trans_when_measurement_is_0 = True 

nonmodeled_proteome_allocation = 0
# nonmodeled_proteome_allocation = 0.503973558
dummy_protein = {'id':'PROSYN-PROTDUMMY','AA abundances':dict()}

path_data = './PinheiroEtAl2020_batchXyl.xlsx'

biom_id = 'BIOSYN-BIOXYL-PINHEIRO'

# set to True to avoid printing .lst files and other large log files
hide_output = True
#doesn't work for GAMS 42 (on Roar cluster): output_redirect_str = ' writeOutput=0 ll=0' if hide_output else ''
#output_redirect_str = ' o=/dev/null ll=0 lo=0' if hide_output else ''
output_redirect_str = ' o=/dev/null' if hide_output else ''

spont_rxn_suffixes = ['SPONT', 'UNKNOWN']
rxns_to_ignore_for_kapps = [] # ignored when printing rxns during C1
