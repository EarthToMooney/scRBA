# %%
# add rxns making products of interest
import pandas as pd
import cobra
from copy import deepcopy

import sys
import os,shutil
sys.path.append('../pycore/')
from utils import extract_details_from_rxnid
from gsm_custom_functions import *

rerun_RBA = True # make False only when testing other parts

# %%
df_rxns = pd.read_excel('./input/pathways.xlsx')
df_prod_info = pd.read_excel('./input/prod_mw.xlsx', header=None)
copy_input_files = False # set to False if you want to customize them
# Manually create input/phenotype_common.txt

# %%
model = cobra.io.load_json_model('../build_model/input/GSM_iSace1144_rba.json')
model = deepcopy(model)
model.solver = 'glpk'
newmodel_path = './input/GSM_iSace1144_rba_pathwayadd.json'

# %%
path_gams = '../GAMS/'
path_app_gams = '../application/input/GAMS_model_application/'
# copy files from GAMS/model folder into path_app_gams, overwriting existing files
if copy_input_files:
	for f in os.listdir(path_gams+'model/'):
		shutil.copy(path_gams+'model/'+f, path_app_gams)

# %%
# model media - same as scRBA
model.objective = dict()

model.reactions.ATPM_c.lower_bound = 1
model.reactions.EX_glc__D_e.bounds = (-13.21,1000)
model.reactions.BIOMASS_AERO_SC_hvd.bounds = (0.1, 0.1)
model.reactions.EX_o2_e.bounds = (-1000,1000)
model.reactions.EX_co2_e.bounds = (-1000,1000)

# Allow YNB media uptake
model.reactions.EX_thm_e.bounds = (-1000,1000)
model.reactions.EX_ribflv_e.bounds = (-1000,1000)
model.reactions.EX_nac_e.bounds = (-1000,1000)
model.reactions.EX_pydxn_e.bounds = (-1000,1000)
model.reactions.EX_fol_e.bounds = (-1000,1000)
model.reactions.EX_pnto__R_e.bounds = (-1000,1000)
model.reactions.EX_4abz_e.bounds = (-1000,1000)
model.reactions.EX_inost_e.bounds = (-1000,1000)

# %%
prod_rxns = dict()
prods = set()

sij_add = set()
added_rxns = set()
new_species = set()

# add rxns to both models
for i in df_rxns.index:
    rxnid = extract_details_from_rxnid(df_rxns.id[i])[1]
    eqn = df_rxns.reaction[i]
    eqn = eqn.replace('MET-', '')
    p = df_rxns.loc[i, 'product']
    
    rxn = cobra.Reaction(rxnid)
    model.add_reactions([rxn])
    rxn.reaction = eqn
    rxn.bounds = (0,0)

    # add to RBA model
    met_dict = metabolites_dict_from_reaction_equation_RBA(rxn.reaction)
    for k,v in met_dict.items():
        if v != '' and k != '':
            sij_add.add("'MET-" + k + "'.'" + df_rxns.id[i] + "' " + str(v))
            new_species.add("'MET-" + k + "'")

    ## Example of how you might count 2 related products together
    if p == '3hpp':
        prod_rxns['3hppa'].append(rxnid)
        prod_rxns['3hppb'].append(rxnid)
    prods.add(p)
    if p in prod_rxns.keys():
        prod_rxns[p].append(rxnid)
    else:
        prod_rxns[p] = [rxnid]
    added_rxns.add(df_rxns.id[i])
prods = list(prods)
print(prods)

# find molecular weights of products
prod_info = dict()
for i in df_prod_info.index:
    if i == 0:
        # add column names to list
        colnames = [df_prod_info.loc[i,j] for j in range(1, len(df_prod_info.loc[i]))]
    else: 
        prod_info[df_prod_info.loc[i,0]] = dict()
        
        # for each column in df_prod_info.loc, add to prod_info with its key being the column title
        for j in range(1, len(df_prod_info.loc[i])):
            prod_info[df_prod_info.loc[i,0]][colnames[j-1]] = df_prod_info.loc[i,j]

# %%
# update sij
def process_GAMS_file(f):
    with open(f, 'r') as f:
        lines = f.readlines()
    # remove all lines that are empty or only contain '/'
    lines = [x.replace('\n','') for x in lines if x.strip() != '/' and x.strip() != '']
    return lines
sij = process_GAMS_file('./input/GAMS_model_application/RBA_sij.txt')
# add new lines
sij += list(sij_add)
# remove duplicates
sij = list(set(sij))
with open('./input/GAMS_model_application/RBA_sij.txt', 'w') as f:
    f.write('\n'.join(['/'] + sij + ['/']))

# add rxns
with open('./input/GAMS_model_application/RBA_rxns_add.txt', 'w') as f:
    f.write('\n'.join(['/'] + list(added_rxns) + ['/']))
allrxns = process_GAMS_file('./input/GAMS_model_application/RBA_rxns.txt')
allrxns += ["'" + x + "'" for x in added_rxns]
allrxns = list(set(allrxns))
with open('./input/GAMS_model_application/RBA_rxns.txt', 'w') as f:
    f.write('\n'.join(['/'] + allrxns + ['/']))

# add prods
species = process_GAMS_file('./input/GAMS_model_application/RBA_species.txt')
# add 'MET-' to products
species += [x for x in new_species]
species = list(set(species))
with open('./input/GAMS_model_application/RBA_species.txt', 'w') as f:
    f.write('\n'.join(['/'] + species + ['/']))

# %%
cobra.io.save_json_model(model, newmodel_path)

# %% [markdown]
#  # Run FBA

# %%
fba_results = {p:dict() for p in prods}

for p in prods:
    if p in ['3hppa', '3hppb','3hpp']:
        exrxn = 'EX_3hpp_e'
        p_for_mw = '3hppa'
    else:
        exrxn = 'EX_' + p + '_e'
        p_for_mw = p
    if p in ['citm', 'polyhb']:
        model.reactions.THRA_c.knock_out()
    for r in prod_rxns[p]:
        model.reactions.get_by_id(r).bounds = (0,1000)
    model.reactions.get_by_id(exrxn).objective_coefficient = 1
    fba = model.optimize()
    
    for r in prod_rxns[p]:
        model.reactions.get_by_id(r).bounds = (0,0)
    if p in ['citm', 'polyhb']:
        model.reactions.THRA_c.bounds = (0,1000)
    model.objective = dict()
    
    fba_results[p]['vprod'] = fba[exrxn]
    fba_results[p]['vglc'] = -fba['EX_glc__D_e']
    fba_results[p]['yield'] = -prod_info[p_for_mw]['MW (g/mol)']*fba[exrxn]/fba['EX_glc__D_e']/180.156 # 180.156 is the molecular weight of glucose
    #print(p, fba.status, fba[exrxn], fba['EX_glc__D_e'])
    print('product exchange flux: '+str(fba[exrxn]))
    print('product exchange flux / glc uptake:',str(fba[exrxn]/fba['EX_glc__D_e']))

# %% [markdown]
#  # Create RBA files from template

# %%
# copy runRBA_max_prod.gms
# if copy_input_files:
# 	for p in prods:
		# copy 

# %% [markdown]
#  # Run RBA

# %%
# shutil.copy(os.path.join(path_gams, 'runRBA_max_prod.gms'),
#             os.path.join(path_out, 'runRBA_max_prod.gms'));
# shutil.copy(os.path.join(path_gams, 'cplex.opt'),
#             os.path.join(path_out, 'cplex.opt'));

if rerun_RBA:
    for p in prods:
            cmds = ['cd output_max/' + p,
                    'module load gams',
                    'python3 runRBA.py']
            os.system('\n'.join(cmds))

# %% [markdown]
#  # Extract RBA yields

# %%
RBA_result_dict = dict()
for prodname in prods:
    if prodname == '3hpp':
        plist = ['3hppa', '3hppb']
    else:
        plist = [prodname]
    for p in plist:
        RBA_result_dict[p] = dict()
        with open('./output_max/' + p + '/report.txt') as f:
            text = f.read().split('\n')

        for t in text[1:]:
            k,v = t.split('\t')
            try:
                RBA_result_dict[p][k] = float(v)
            except:
                RBA_result_dict[p][k] = v
            
        print(RBA_result_dict[p]['vprod'])

# %%
# copy compiled_results_template.xlsx to compiled_results.xlsx
shutil.copy(os.path.join('./compiled_results_template.xlsx'),
			os.path.join('./compiled_results.xlsx'))
# write FBA yields to rateFBA column (2nd how has headers)
df = pd.read_excel('./compiled_results.xlsx', header=1)
for prodname in prods:
    if prodname == '3hpp':
        plist = ['3hppa', '3hppb']
    else:
        plist = [prodname]
    for p in plist:
        # add each column from prod_info to df
        for k,v in prod_info[p].items():
            df.loc[df['product'] == p, k] = v
        df.loc[df['product'] == p, 'rateFBA'] = fba_results[p]['vprod']
        df.loc[df['product'] == p, 'vglcFBA'] = fba_results[p]['vglc']
        df.loc[df['product'] == p, 'yieldFBA'] = fba_results[p]['yield']
        df.loc[df['product'] == p, 'rateRBA'] = RBA_result_dict[p]['vprod']
        df.loc[df['product'] == p, 'vglcRBA'] = RBA_result_dict[p]['vglc']
        df.loc[df['product'] == p, 'yieldRBA'] = RBA_result_dict[p]['yield']
        df.loc[df['product'] == p, 'yRBA/yFBA'] = RBA_result_dict[p]['yield']/fba_results[p]['yield']*100

df.to_excel('./compiled_results.xlsx', index=False)
