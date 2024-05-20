# %%
import pandas as pd
import sys
sys.path.append('../../../../pycore/')
from simulate import RBA_result
from utils import extract_details_from_rxnid

# %%
# Load enzyme info
df_enz = pd.read_excel('../../../../build_model/input/ENZYME_stoich_curation_2021-10-01.xlsx')

# Load path of set4 (enzyme-reaction many-to-many mapping)
set4_path = '../kapp_ambiguousLoad_case_resolve_common.txt'

# %%
biom_id = 'BIOSYN-BIODILNOGAM'

res_metab = RBA_result(biom_id=biom_id)
res_metab.load_raw_flux('./min_flux_sum/min_flux_sum.flux.txt')
res_metab.calculate_metabolic_flux()

res_esyn = RBA_result(biom_id=biom_id, twocol_format=True)
res_esyn.load_raw_flux('./enz_from_proteome/enz_flux_calculation.txt')

mu = res_metab.growth_rate
print('Growth rate:', mu)

# %%
print('EX_glc__D_e', res_metab.metabolic_flux['EX_glc__D_e'])
print('EX_o2_e', res_metab.metabolic_flux['EX_o2_e'])
for k,v in res_metab.metabolic_flux.items():
    if k[:3] == 'EX_' and v > 1e-6:
        print(k,v)

# %%
#### Map rxn to enz
rxndict = {k:[] for k,v in res_metab.metabolic_flux.items() if abs(v) > 0}

enzdict = dict()
for k,v in res_esyn.raw_flux.items():
    if k.split('-')[0] == 'ENZSYN' and v > 0:
        _,enz = k.split('-', maxsplit=1)
        enzdict[enz] = []
        
for i in df_enz.index:
    rxn = df_enz.rxn_src[i]
    enz = df_enz.enz[i]
    if rxn in rxndict.keys() and enz in ['SPONT', 'UNKNOWN']:
        rxndict[rxn].append('zeroCost')
        
    if rxn in rxndict.keys() and enz in enzdict.keys():
        rxndict[rxn].append(enz)
        enzdict[enz].append(rxn)
        
rxndict = {k:set(v) for k,v in rxndict.items()}
rxndict = {k:v for k,v in rxndict.items() if v != {'zeroCost'}}
enzdict = {k:set(v) for k,v in enzdict.items()}

# %%
# Find enzymes whose together carry a total load of reactions
x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
enz_share_rxn_load = set().union(*[v for v in x.values()])

# Find enzymes whose individually carry loads of multiple reactions
x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
enz_multiload = set(x.keys())

# Set 1: Enzyme-reaction one-to-one load mapping
set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])

# Set 2: Enzyme-reaction one-to-many
set2 = enz_multiload - enz_share_rxn_load

# Set 3: Enzyme-reaction many-to-one
set3 = enz_share_rxn_load - enz_multiload

# Set 4: Enzyme-reaction many-to-many
set4 = enz_share_rxn_load & enz_multiload

# %%
# Check manual resolve of ambiguous load case covers set4
with open(set4_path) as f:
    text = f.read().split('\n')[1:]
text = [i for i in text if i != '']
enzs = []
for i in text:
    enzs += i.split('\t')[0].split(',')
for i in set4 - set(enzs):
    # Manually verified these enzymes catalyzing very low flux values / unverified flux
    # Ignore in kapp calculation
    if i not in ['YOR348C', 'YMR313C']:
        print(i)

# %%
# Find enzymes whose together carry a total load of reactions
x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
enz_share_rxn_load = set().union(*[v for v in x.values()])

# Find enzymes whose individually carry loads of multiple reactions
x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
enz_multiload = set(x.keys())

# Set 1: Enzyme-reaction one-to-one load mapping
set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])

# Set 2: Enzyme-reaction one-to-many
set2 = enz_multiload - enz_share_rxn_load

# Set 3: Enzyme-reaction many-to-one
set3 = enz_share_rxn_load - enz_multiload

# Set 4: Enzyme-reaction many-to-many
with open(set4_path) as f:
    text = f.read().split('\n')[1:]
text = [i for i in text if i != '']
enzs = []
for i in text:
    enzs += i.split('\t')[0].split(',')
set4 = set(enzs)
    
# In manual checking, some enzymes in set2 and set3 requires special treatment in calculation
# which are recorded in set4. Thus excluded those in set2 and set3
set2 = set2 - set4
set3 = set3 - set4

kapp = dict()

### Set 1: Enzyme-reaction one-to-one load mapping
for enz in set1:
    rxn = [i for i in enzdict[enz]][0]
    rval = res_metab.metabolic_flux[rxn]
    if rval > 0:
        rdir = 'FWD'
    elif rval < 0:
        rdir = 'REV'
    else:
        print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)
    
    rid = 'RXN-' + rxn + '_' + rdir + '-' + enz
    
    enzval = res_esyn.raw_flux['ENZSYN-' + enz]
    kapp[rid] = mu * abs(rval) / enzval / 3600
    
### Set 2: Enzyme-reaction one-to-many
for enz in set2:
    rids = []; rvalsum = 0;
    for rxn in enzdict[enz]:
        rval = res_metab.metabolic_flux[rxn]
        if rval > 0:
            rdir = 'FWD'
        elif rval < 0:
            rdir = 'REV'
        else:
            print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)

        rids.append('RXN-' + rxn + '_' + rdir + '-' + enz)
        rvalsum += abs(rval)
        
    enzval = res_esyn.raw_flux['ENZSYN-' + enz]
    
    for rid in rids:
        kapp[rid] = mu * rvalsum / enzval / 3600
        
### Set 3: Enzyme-reaction many-to-one
rxns = set().union(*[enzdict[enz] for enz in set3])
for rxn in rxns:
    rval = res_metab.metabolic_flux[rxn]
    if rval > 0:
        rdir = 'FWD'
    elif rval < 0:
        rdir = 'REV'
    else:
        print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)
        
    rids = []; enzvals = [];
    for enz in rxndict[rxn]:
        rids.append('RXN-' + rxn + '_' + rdir + '-' + enz)
        enzvals.append(res_esyn.raw_flux['ENZSYN-' + enz])
    enzval = max(enzvals)
        
    for rid in rids:
        kapp[rid] = mu * abs(rval) / enzval / 3600
        
### Set 4: Enzyme-reaction many-to-many
with open(set4_path) as f:
    text = f.read().split('\n')[1:]
text = [i for i in text if i != '']
mapper = dict()
for line in text:
    k,v,rule = line.split('\t')
    if v not in ['noRxnFlux_resolveByManualCheck', 'disableManually']:
        mapper[k] = (v,rule)
        
for enztext,x in mapper.items():
    rxntext,rule = x
    
    # Parsing enzymes
    enzs = enztext.split(',')
    if rule == 'divmaxenz':
        enzs_measured = []
        for enz in enzs:
            if 'ENZSYN-' + enz in res_esyn.raw_flux.keys():
                enzs_measured.append(enz)
        if enzs_measured == []:
            continue
        else:
            enzval = max([res_esyn.raw_flux['ENZSYN-'+enz] for enz in enzs_measured])
            if enzval == 0:
                continue
            
    else:
        enzval = 0
        enzs_measured = []
        for enz in enzs:
            if 'ENZSYN-' + enz in res_esyn.raw_flux.keys():
                enzs_measured.append(enz)
                enzval += res_esyn.raw_flux['ENZSYN-'+enz]
        if enzval == 0:
            continue
    
    # Parsing reactions
    rxns = rxntext.split(',')
    rxnval = 0
    rxns_on = []
    for rxn in rxns:
        if rxn in res_metab.metabolic_flux.keys():
            rxns_on.append(rxn)
            rxnval += abs(res_metab.metabolic_flux[rxn])
    if rxnval < 1e-5:
        continue
    
    rids = []
    for rxn in rxns_on:
        rval = res_metab.metabolic_flux[rxn]
        if rval > 0:
            rdir = 'FWD'
        elif rval < 0:
            rdir = 'REV'
        else:
            print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)
            
        for enz in enzs_measured:
            rids.append('RXN-' + rxn + '_' + rdir + '-' + enz)
            
    for rid in rids:
        kapp[rid] = mu * rxnval / enzval / 3600

# %% [markdown]
# #### Write calculated kapp (1/s)

# %%
#Flux is numerically low, near zero
with open('../exclude_parameterization_list.txt') as f:
    excl = f.read().split('\n')
excl = [r for r in excl if r != '']
excl += []

texts = ['rxnid\tkapp (1/s)']
perhr_texts = ['rxnid\tkapp (1/s)']
for k,v in kapp.items():
    if k not in excl:
        texts.append(k + '\t' + str(v))
        perhr_texts.append(k + '\t' + str(v*3600))
    
with open('./kapps_in_vivo.txt', 'w') as f:
    f.write('\n'.join(texts))
with open('./kapps_per_hr.txt', 'w') as f:
    f.write('\n'.join(texts))

# %%
import numpy as np
np.median(list(kapp.values()))

# %%
