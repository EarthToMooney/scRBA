import pandas as pd
import sys
sys.path.append('../../../../pycore/')
from simulate import RBA_result
from utils import extract_details_from_rxnid
import os
# Load enzyme info
df_enz = pd.read_excel('../../../../build_model/input/ENZYME_stoich_curation_2021-10-01.xlsx')

# Load path of set4 (enzyme-reaction many-to-many mapping)
set4_path = '../kapp_ambiguousLoad_case_resolve_common.txt'
# Eric: How this file works: handles cases of multiple enzymes corresponding to multiple rxns (many-to-many)
# columns: 
#   enzyme (for all enz that correspond to a rxn)
#   rxn: 
#       disableManually for any I don't want to count
#       noRxnFlux_resolveByManualCheck: behaves identically to disableManually, but perhaps serves as a warning
#   rule: 
#       divmaxenz for cases where I want to treat the [enzyme] as equaling its maximum value
biom_id = 'BIOSYN-BIODILNOGAM'

res_metab = RBA_result(biom_id=biom_id)
res_metab.load_raw_flux('./min_flux_violation/enz_alloc.flux.txt')
#res_metab.load_raw_flux('./min_flux_violation/min_flux_violation.flux.txt')
res_metab.calculate_metabolic_flux()

res_esyn = RBA_result(biom_id=biom_id, twocol_format=True)
res_esyn.load_raw_flux('./min_flux_violation/enz_flux_calculation.txt')

with open("../../../../GAMS/model/RBA_rxns_rxnmetabolicnetwork.txt","r") as f: 
    all_rxns = [line.strip().replace("'","") for line in f.readlines() if line.strip() != "/"]

mu = res_metab.growth_rate
print('Growth rate:', mu)
print('EX_glc__D_e', res_metab.metabolic_flux['EX_glc__D_e'])
print('EX_o2_e', res_metab.metabolic_flux['EX_o2_e'])
for k,v in res_metab.metabolic_flux.items():
    if k[:3] == 'EX_' and v > 1e-6:
        print(k,v)
#### Map rxn to enz
rxndict = {k:[] for k,v in res_metab.metabolic_flux.items() if abs(v) > 0}

enzdict = dict() # all enzymes made by the model
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
rxndict_zeroCost = {k:v for k,v in rxndict.items() if v == {'zeroCost'}}
rxndict = {k:v for k,v in rxndict.items() if v != {'zeroCost'}}
enzdict = {k:set(v) for k,v in enzdict.items()}
# Find enzymes that together carry a total load of reactions
x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
enz_share_rxn_load = set().union(*[v for v in x.values()])

# Find enzymes that individually carry loads of multiple reactions
x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
enz_multiload = set(x.keys())

# Set 1: Enzyme-reaction one-to-one load mapping
set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])
# write set1 to file
with open('./kapp_set1_one_to_one.txt', 'w') as f:
    f.write('\n'.join(set1))

# Set 2: Enzyme-reaction one-to-many
set2 = enz_multiload - enz_share_rxn_load
# write set2 to file
with open('./kapp_set2_one_to_many.txt', 'w') as f:
    f.write('\n'.join(set2))

# Set 3: Enzyme-reaction many-to-one
set3 = enz_share_rxn_load - enz_multiload
# write set3 to file
with open('./kapp_set3_many_to_one.txt', 'w') as f:
    f.write('\n'.join(set3))

# Set 4: Enzyme-reaction many-to-many
set4 = enz_share_rxn_load & enz_multiload
# write set4 to file
with open('./kapp_set4_many_to_many.txt', 'w') as f:
    f.write('\n'.join(set4))
# if set4_path doesn't exist or is only 1 line long, create it with the header
if not os.path.exists(set4_path) or len(open(set4_path).readlines()) < 2:
    with open(set4_path, 'w') as f:
        f.write('Enzyme\tReaction\tRule\n')
        # add enzymes from set4 to set4_path, and all rxns they catalyze
        for enz in set4:
            f.write(enz + '\t' + ','.join(list(enzdict[enz])) + '\t\n')

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
# Find enzymes whose together carry a total load of reactions
x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
enz_share_rxn_load = set().union(*[v for v in x.values()])

# Find enzymes whose individually carry loads of multiple reactions
x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
enz_multiload = set(x.keys())

# Set 1: Enzyme-reaction one-to-one load mapping
set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])

# Set 2: Enzyme-reaction one-to-many *(see note below)
set2 = enz_multiload - enz_share_rxn_load

# Set 3: Enzyme-reaction many-to-one *(see note below)
set3 = enz_share_rxn_load - enz_multiload

# Set 4: Enzyme-reaction many-to-many
with open(set4_path) as f:
    text = f.read().split('\n')[1:]
text = [i for i in text if i != '']
enzs = []
for i in text:
    enzs += i.split('\t')[0].split(',')
set4 = set(enzs)
    
# From manual checking of the sets, some enzymes in set2 and set3 requires special treatment in calculation
# which are recorded in set4. Thus excluded those in set2 and set3
set2 = set2 - set4
set3 = set3 - set4

kapp = dict()
kapp_per_hr = dict()
### For testing only: set each rxn-enzyme's kapp individually
kapp_minimal_assumptions = dict()
kapp_ma_cutoff = 1e-1
kapp_ma_text = []
# print(enzdict)
for enz in enzdict:
    for rxn in enzdict[enz]:
# for enz,rxn in enzdict.items():
    # print(enz,rxn)
    # if rxn isn't an empty set
        if rxn != set():
            if rxn in rxndict.keys():
                rval = res_metab.metabolic_flux[rxn]
                if rval > 0:
                    rdir = 'FWD'
                elif rval < 0:
                    rdir = 'REV'
                else:
                    print('rval == 0, check enzyme ' + enz + ' and reaction ' + rxn)
                
                rid = 'RXN-' + rxn + '_' + rdir + '-' + enz
                # check if raw flux is available
                if 'ENZLOAD-' + rxn + '_' + rdir + '-' + enz in res_esyn.raw_flux.keys():
                    enzval = res_esyn.raw_flux['ENZLOAD-' + rxn + '_' + rdir + '-' + enz]
                    if enzval > 0:
                        kapp_minimal_assumptions[rid] = mu * abs(rval) / enzval
                        if kapp_minimal_assumptions[rid] > kapp_ma_cutoff: # cutoff to prevent low kapps from being used, in case they cause issues
                            kapp_ma_text.append("'" + rid + "'\t" + str(kapp_minimal_assumptions[rid]))

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
    kapp_per_hr[rid] = kapp[rid] * 3600
    
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
        kapp_per_hr[rid] = kapp[rid] * 3600
        
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
    enzval = sum(enzvals) # changed to match scRBA suppMat description
        
    for rid in rids:
        kapp[rid] = mu * abs(rval) / enzval / 3600
        kapp_per_hr[rid] = kapp[rid] * 3600
        
        
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
            # changed to match scRBA suppMat description (pg. 16 of suppMat1)
            enzval = sum([res_esyn.raw_flux['ENZSYN-'+enz] for enz in enzs_measured]) 
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
        if rxn not in all_rxns:
            print('Reaction from '+set4_path+' ignored due to its absence from model: '+rxn)
            continue
        if rxn in res_metab.metabolic_flux.keys():
            rxns_on.append(rxn)
            rxnval += abs(res_metab.metabolic_flux[rxn])
    if rxnval < 1e-5:
        print('total flux below cutoff for ' + enztext + ' and ' + rxntext)
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
        kapp_per_hr[rid] = kapp[rid] * 3600
#Flux is numerically low, near zero
with open('../exclude_parameterization_list.txt') as f:
    excl = f.read().split('\n')
# # TEST: add all rxns with fluxes below 1e-5 to excl
# for k,v in res_metab.metabolic_flux.items():
#    if abs(v) < 1e-5:
#        excl.append(k)
excl = [r for r in excl if r != '']
excl += []

import numpy as np
default_kapp = np.median(list(kapp.values())) * 3600

texts = ['rxnid\tkapp (1/s)']
perhr_texts_used_only = ['/']
kapp_per_hr = dict()
for k,v in kapp.items():
    if k not in excl:
        texts.append(k + '\t' + str(v))
        kapp_per_hr[k] = v*3600
    
with open('./kapps_in_vivo.txt', 'w') as f:
    f.write('\n'.join(texts))
kapp_max = max(list(kapp.values())) * 3600
kapp_ma_med = np.median(list(kapp_minimal_assumptions.values()))
# return kapps in a format suitable for use in RBA
kapp_txt = ['/']

# check if any inactive rxns were used in last step
with open('min_flux_violation/min_flux_violation.flux_essential_inactive_rxns.txt') as f:
    rxns_essential_inactive = f.read().split('\n')
rxns_essential_inactive = [i for i in rxns_essential_inactive if i != '']
rxns_essential_inactive = [i.split('\t')[0] for i in rxns_essential_inactive]

output_info = []
for r in rxns_essential_inactive:
    if r not in kapp_per_hr.keys():
        kapp_per_hr[r] = default_kapp
for k,v in kapp_per_hr.items():
    perhr_texts_used_only.append("'" + k + "'" + '\t' + str(v))

with open('./kapps_per_hr_without_unused_rxns.txt', 'w') as f:
    f.write('\n'.join(perhr_texts_used_only + ['/']))

perhr_texts = ['/']
# add all other rxns with default kapp
for rxn in all_rxns: 
    enz = rxn.split('-')[-1]
    rxn_name_without_enz = rxn.split('-')[-2]
    if enz not in ['SPONT', 'UNKNOWN']:
        if rxn in kapp_per_hr.keys():
            pass
            # if kapp_per_hr[rxn] <= 1e-5:
                # output_info.append('kapp <= cutoff for ' + rxn)
                # kapp_per_hr[rxn] = default_kapp
        elif rxn in rxns_essential_inactive:
            kapp_per_hr[rxn] = default_kapp # to minimize the risk of kapps making growth infeasible 
        else: 
            output_info.append('No kapp found for ' + rxn)
            kapp_per_hr[rxn] = default_kapp
        if rxn not in kapp_minimal_assumptions.keys():
            # set kapp to the max value of all kapps for that enzyme (or reaction, if unavailable), to minimize risk of overestimation
            new_kapp = 0
            for k,v in kapp_minimal_assumptions.items():
                if enz in k.split('-')[-1]:
                    if v > new_kapp:
                        new_kapp = v
                # find max kapp for that reaction
                elif rxn_name_without_enz == k.split('-')[-2]:
                    if v > new_kapp:
                        new_kapp = v
                # remove location to see if enzyme is identical to another
                elif enz == k.split('-')[-2].split('_')[-1]:
                    if v > new_kapp:
                        new_kapp = v
            if new_kapp == 0:
                output_info.append('No kapp found for ' + rxn + ' or enzyme ' + enz)
                new_kapp = kapp_ma_med
            kapp_minimal_assumptions[rxn] = new_kapp
            kapp_ma_text.append("'" + rxn + "'\t" + str(kapp_minimal_assumptions[rxn]))
        elif kapp_minimal_assumptions[rxn] <= kapp_ma_cutoff:
            kapp_ma_text.append("'" + rxn + "'\t" + str(kapp_ma_med))
        kapp_txt.append("'" + rxn + "'\t" + str(kapp_per_hr[rxn]))
        perhr_texts.append("'" + rxn + "'\t" + str(kapp_per_hr[rxn]))

with open('./kapps_per_hr.txt', 'w') as f:
    f.write('\n'.join(perhr_texts + ['/']))
with open('./kapp_default_values.txt', 'w') as f:
    f.write('\n'.join(output_info))
with open('./kapps_RBA.txt', 'w') as output:
    output.write('\n'.join(kapp_txt) + '\n' + '/')
with open('./kapps_minimal_assumptions.txt', 'w') as output:
    output.write('/' + '\n' + '\n'.join(kapp_ma_text) + '\n' + '/')
# if any rxns are inactive, we've predicted their kapps and know their fluxes
# use these to infer ENZSYN fluxes (ENZSYN=mu*flux/kapp)
kapp_test_text = []
if len(rxns_essential_inactive) > 0:
    # # load fluxes
    # with open('min_flux_violation/min_flux_violation.flux.txt') as f:
    #     fluxes = f.read().split('\n')
    # fluxes = [i.split('\t') for i in fluxes]
    # fluxes = {i[0]:float(i[1]) for i in fluxes}
    # # calculate ENZSYN fluxes
    # enzsyn_fluxes = []
    c = 0
    for rxn in rxns_essential_inactive:
        if rxn in kapp_per_hr.keys():
            # assume kapp is at its max value, to prevent excessive constraints
            kapp_per_hr[rxn] = kapp_max 
            kapp_test_text.append('Equation EnzCap'+str(c)+'; EnzCap'+str(c)+".. v('"+rxn.replace('RXN-','ENZLOAD-')+"')*"+str(kapp_per_hr[rxn])+" =e= %mu% * v('"+rxn+"');")
            c += 1
        # add its kapp as a constraint on ENZLOAD and flux
    #     if rxn in kapp_per_hr.keys() and rxn in fluxes.keys():
    #         enzsyn_fluxes.append('ENZSYN-' + rxn + '\t' + str(mu*fluxes[rxn]/kapp_per_hr[rxn]))
    # with open('./enzsyn_fluxes.txt', 'w') as f:
    #     f.write('\n'.join(enzsyn_fluxes))
    with open('./kapp_test.txt', 'w') as f:
        f.write('\n'.join(kapp_test_text))
