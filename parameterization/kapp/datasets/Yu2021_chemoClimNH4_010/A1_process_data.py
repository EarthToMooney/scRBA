# %%
import pandas as pd
import numpy as np

import sys
sys.path.append('../../../../pycore/')
from utils import metabolites_dict_from_reaction_equation_RBA

# %%
mu = 0.1
cols_data = ['prot.4', 'prot.5', 'prot.6']
raw_data_path = '../raw_data_files/Yu2021_rawdata.xlsx'
fout_name = './Yu2021_chemoClimNH4_010.xlsx'

# %%
# Load path
path_gen = '../../../../build_model/'

prot_path = path_gen + 'input/PROTEIN_stoich_curation_2021-09-28.xlsx'
model_xlsx_path = path_gen + 'model/RBA_stoichiometry.xlsx'
ribonuc_path = path_gen + 'input/RIBOSOME_nucleus.xlsx'
ribomito_path = path_gen + 'input/RIBOSOME_mitochondria.xlsx'

# %%
verified_dubious_genes = ['YJL182C', 'YPL044C', 'YML009W-B']
df_raw = pd.read_excel(raw_data_path, sheet_name='Supplementary File 1b', skiprows=[0])
df_raw.index = df_raw['standard_id'].to_list()
idx = [i for i in df_raw.index if i not in verified_dubious_genes]
df_raw = df_raw.loc[idx, :]

# Load protein
df_prot = pd.read_excel(prot_path)
df_prot.index = df_prot.id.to_list()
# Strip compartment
df_prot.index = [i.split('_')[0] if '_' in i else i for i in df_prot.index]
df_prot['id'] = df_prot.index.to_list()
df_prot = df_prot.drop_duplicates(subset=['id'])
# Protein copy selector: start with empty file
# per Hoang: "protein_copies_selector.txt" is a manually written file. It should be empty for you in the first run. Then, you have to check your kapp calculation to see if you need to add any entries. For example, an enzyme ILV2-ILV6 complex and ILV2 are functional, but I choose to assign estimated kapp value to ILV2-6 complex only (which I think is the primary complex for catalysis, otherwise it doesn't make sense to me why ILV2-6 complex is even needed). There will be many cases like this where you have to write entries to the "protein_copies_selector.txt" file
df_select = pd.read_csv('./input/protein_copies_selector.txt', sep='\t')
df_select.index = df_select.gene_src.to_list()

# Ribosome (nucleus and mitochondrial)
df_ribonuc = pd.read_excel(ribonuc_path)
df_ribomito = pd.read_excel(ribomito_path)

# %%
#### HANDLE MISSING MEASUREMENTS FOR SUBUNIT COMPONENT OF HETEROMERIC ENZYMES
# E.g., missing subunit measurements for ATP synthase complex
# Stoichiometry
df_eqn = pd.read_excel(model_xlsx_path)
df_eqn.index = df_eqn.id.to_list()

# %%
#### MW
df_mw = pd.read_csv('../scProteins_MW.csv', sep='\t')
df_mw.index = df_mw.gene_id.to_list()

# %% [markdown]
# Process data

# %%
pdata_raw = df_raw.loc[:, cols_data].mean(axis=1).fillna(0)
mw = df_mw.loc[pdata_raw.index.to_list(), 'MW (g/mmol)']

weight_tot = sum(pdata_raw * mw)
# pdata variable: protein fraction in proteome (g protein / gDW)
pdata = (pdata_raw * mw) / weight_tot
pdata = pdata[pdata > 0]

# %%
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)', 'type', 'conc (g/gDW)', 'vtrans (mmol/gDW/h)']
idx = [i for i in df_prot.index if i in pdata.index]

df_data = pd.DataFrame(columns=cols, index=idx)
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)']
df_data.loc[idx, cols] = df_prot.loc[idx, cols]

ptot = (36.94 + 34.22*mu) / 100 #Clim
#ptot = (10.57 + 108.56*mu) / 100 #Nlim

for i in df_data.index:
    pval = pdata[i]
    mw = df_data.loc[i, 'MW (g/mmol)']
    df_data.loc[i, 'conc (g/gDW)'] = pval * ptot
    df_data.loc[i, 'vtrans (mmol/gDW/h)'] = mu * pval * ptot / mw
    df_data.loc[i, 'type'] = 'truedata_enz'
    
    if i in df_ribonuc.id.to_list():
        df_data.loc[i, 'type'] = 'truedata_ribonuc'
    elif i in  df_ribomito.id.to_list():
        df_data.loc[i, 'type'] = 'truedata_ribomito'

# %%
# Store index that matches the raw data
idx_truedata_old = df_data.index.to_list()

# Reindex - incorporate info from protein copy selector
idx = [df_select.selected_compartmental_copy[i] if i in df_select.index \
       else i for i in df_data.index]
df_data.index = idx
df_data['id'] = df_data.index.to_list()

# %%
# Clean out NaN rows
df_data = df_data[df_data['conc (g/gDW)'].isnull() == False]

# %% [markdown]
# Gap-fill data

# %%
# Load protein
df_prot = pd.read_excel(prot_path)
df_prot.index = df_prot.id.to_list()

# %%
idx_enzsyn = df_eqn[df_eqn.id.str.contains('ENZSYN-')].index
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)']

for i in idx_enzsyn:
    x = metabolites_dict_from_reaction_equation_RBA(df_eqn.reaction[i])
    met_dict = dict()
    for k,v in x.items():
        if k == '':
            continue
        if v.is_integer():
            met_dict[k] = int(v)
        else:
            met_dict[k] = v

# %%
idx_enzsyn = df_eqn[df_eqn.id.str.contains('ENZSYN-')].index
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)']

for i in idx_enzsyn:
    x = metabolites_dict_from_reaction_equation_RBA(df_eqn.reaction[i])
    met_dict = dict()
    for k,v in x.items():
        if k == '':
            continue
        if v.is_integer():
            met_dict[k] = int(v)
        else:
            met_dict[k] = v
            
    met_dict = {k.split('-', maxsplit=1)[1]:v for k,v in met_dict.items() if v < -1e-6}
    in_data = set(met_dict) & set(df_data.index)
    if len(in_data) > 0.5 and len(in_data) < len(met_dict):
        #print(i, len(in_data), len(met_dict), ','.join(in_data))
        vmin = min([df_data.loc[k, 'vtrans (mmol/gDW/h)'] / met_dict[k] for k in met_dict.keys() if k in in_data])
        cmin = min([df_data.loc[k, 'conc (g/gDW)'] / met_dict[k] for k in met_dict.keys() if k in in_data])
        for k in met_dict.keys():
            if k not in in_data:
                idx_truedata_old.append(k)
                df_data.loc[k, cols] = df_prot.loc[k, cols]
                df_data.loc[k, 'conc (g/gDW)'] = cmin * met_dict[k]
                df_data.loc[k, 'vtrans (mmol/gDW/h)'] = vmin * met_dict[k]
                df_data.loc[k, 'type'] = 'gapfill_subunit'

# %%
df_data.to_excel(fout_name, index=None)

# %% [markdown]
# Calculate fraction of non-enzymatic and non-ribosomal proteins

# %%
idx = [i for i in pdata.index if i not in idx_truedata_old]
pdata[idx].sum()
