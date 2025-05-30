# Model-specific settings for kapp calculations, as well as steps that must be run each time
import pandas as pd
import numpy as np
import requests,json,os,shutil,sys
from pycore.cobrapy_custom_extras import *

model_root_path = '../../../../'
path_model = model_root_path+'GAMS/model/'
pycore_path = model_root_path+'pycore/'
for path in [model_root_path, pycore_path]:
    if path not in sys.path:
        sys.path.append(path)
# retrieve name of the directory containing the options file
dir_name = os.path.basename(os.path.dirname(os.path.abspath(__file__)))

from pycore.utils import metabolites_dict_from_reaction_equation_RBA, extract_details_from_rxnid
from pycore.simulate import * 
from pycore.gsm_custom_functions import *

vmax = 1e3 # max flux in either direction, before applying nscale
nscale=1e5 # will be used as scale when running GAMS
nscale2=1e5 # 2nd scale, if different one needed; otherwise, set to same as nscale

# Load path
path_gen = model_root_path+'build_model/'
path_gams = model_root_path+'GAMS/'

sm_path = path_gen + 'input/GSM_iSace1144_rba.json'

prot_path = path_gen + 'model/PROTEIN_stoich_curation.tsv'
unknown_prot_path = path_gen + 'input/PROTEIN_stoich_curation_unknown.xlsx'
model_xlsx_path = path_model + 'RBA_stoichiometry.tsv'
metab_rxns_path = path_model + 'RBA_rxns_rxnmetabolicnetwork.txt'
gsm_rxn_ids_path = path_model + 'SM_rxn_ids.txt'
sij_path = path_model + 'RBA_sij.txt'
aa_mapping_path = path_gen + 'input/PROTEIN_amino_acid_map.txt'
aa_map = pd.read_csv(aa_mapping_path, sep='\t')
aa_dict = dict(zip(aa_map['aa_abbv'], aa_map['MW']))

# for flux data; keep path blank if none provided
flux_path = ''
col_LB = 'mfaLB'
col_UB = 'mfaUB'

# must match the growth rate in other files 
# 	may help to add or subtract measurement error in the value itself, for future reference (e.g., mu = 0.060 - 0.001, not mu = 0.059)
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
proteomics_units = 'g/g protein' # supported options: 'g/g protein' (originally built with this in mind), 'g/gDW'
use_ribo_data = False # False to ignore ribosome subunit abundance data # False to ignore ribosome subunit abundance data (still used for MPFS constraints)
find_kribo = True # if true, runs kribo bisection

# if True, proteins with measured translation rates of 0 can be made at a flux equal to prosynSlackUB, thus discouraging their use where possible. If you want such proteins to be made, you thus may have to specify them in the phenotype.txt file as having fixed values for prosynSlackUB. 
# if False, proteins with measured translation rates of 0 can't be made. This may cause issues if no isozymes exist for such proteins but their rxns are essential, and such proteins may still be present but unmeasured due to measurement errors.
allow_trans_when_measurement_is_0 = True 

nonmodeled_proteome_allocation = 0
# nonmodeled_proteome_allocation = 0.503973558
dummy_protein = {'id':'PROSYN-PROTDUMMY','AA abundances':dict()}

path_data = './protein_data.tsv'
if os.path.exists(path_data):
    df_data = read_spreadsheet(path_data)
    df_data.index = df_data['id'].to_list()
    df_data = df_data[df_data['conc (g/gDW)'] > 0]
    # Excluding ribosome protein subunit (conflicting if fit to both enzymatic and ribosomal protein data)
    # make perfect copy of df_data
    df_data_full = df_data.copy()
    if not use_ribo_data:
        df_data = df_data[(df_data.type == 'truedata_enz') | (df_data.type == 'gapfill_subunit')]

biom_id = 'BIOSYN-BIOXYL-PINHEIRO'

# set to True to avoid printing .lst files and other large log files
hide_output = True
#doesn't work for GAMS 42 (on Roar cluster): output_redirect_str = ' writeOutput=0 ll=0' if hide_output else ''
#output_redirect_str = ' o=/dev/null ll=0 lo=0' if hide_output else ''
output_redirect_str = ' o=/dev/null' if hide_output else ''

spont_rxn_suffixes = ['SPONT', 'UNKNOWN']
rxns_to_ignore_for_kapps = [] # ignored when printing rxns during C1
