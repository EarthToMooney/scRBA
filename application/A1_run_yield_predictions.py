# add rxns making products of interest
import pandas as pd
import cobra
from copy import deepcopy
import json
from cobra.io import load_model

import sys,os,shutil
sys.path.append('../pycore/')
from gsm_custom_functions import *

# (optional) ensures output is visible when running as slurm job
from functools import partial
print = partial(print, flush=True)

# Unless you force all uptake fluxes to specific values, or somehow make yields into an objective function or constraint (i.e., linear fractional programming), they may vary.
# The minYield constraint addresses this, so keep it in.
rerun_FBA = True # make False only when testing other parts.
rerun_RBA = True # make False only when testing other parts. Don't need to turn off on your system if you can't run it, since it'll automatically avoid running RBA in that case.
find_min_media_yields = False # turned off since I use RBA for minimal media where needed
adjust_constraints_if_infeas = False # allow adjustment of constraints (e.g., growth rate, substrate uptake) if production is infeasible.
copy_input_files = True # set to False if you want to customize them
# Manually create input/phenotype_common.txt

df_rxns = pd.read_excel('./input/pathways.xlsx')
df_prod_info = pd.read_excel('./input/prod_mw.xlsx', header=None)

modelname = 'GSM_iSace1144_rba'
model = cobra.io.load_json_model('../build_model/input/'+modelname+'.json')
model = deepcopy(model)
model.solver = 'glpk'
no_gam = False # WIP: add way for FBA to be affected by this
RBA_biomass_rxn = 'BIOSYN-BIODILAERO-NOGAM' if no_gam else 'BIOSYN-BIODILAERO'
FBA_biomass_rxn = 'BIOMASS_MFA_NO_GAM' if no_gam else 'BIOMASS_AERO_SC_hvd'
newmodel_path = './input/'+modelname+'_pathwayadd.json'
vmax = 1000 # highest absolute value of any flux in GSM model

outputFolder = './output_max_withoutMP/'
try:
    # fba_results,RBA_result_dict = dict(),dict()
    fba_results = json.load(open(fba_results_path))
    RBA_result_dict = json.load(open(rba_results_path))
except:
    pass

path_gams = '../GAMS/'
path_app_gams = '../application/input/GAMS_model_application/'
# copy files from GAMS/model folder into path_app_gams, overwriting existing files
if copy_input_files:
	for f in os.listdir(path_gams+'model/'):
		shutil.copy(path_gams+'model/'+f, path_app_gams)

prod_rxns = dict()
prods = set()

sij_add = set()
added_rxns = set()
added_rxns_dict = dict()
new_species = set()

# if run with command line arguments, use them to update variables
if len(sys.argv) > 1:
    # update each variable in the argument list
    for i in range(1,len(sys.argv)):
        # split the argument into key and value
        k,v = sys.argv[i].split('=') # e.g., 'rerun_FBA=True'
        # if the value is 'True' or 'False', convert it to a boolean
        if v == 'True':
            v = True
        elif v == 'False':
            v = False
        # update the variable
        exec(k + ' = ' + str(v))
        print('variable',k,'set to',str(v))
fba_results_path = outputFolder+'all_fba_results.json'
rba_results_path = outputFolder+'all_RBA_results.json'

# add rxns to both models
for i in df_rxns.index:
    tag,rxnid,rxn_dir,enz_id = extract_details_from_rxnid(df_rxns.id[i])
    eqn = df_rxns.reaction[i]
    eqn = eqn.replace('MET-', '')
    p = df_rxns.loc[i, 'product']
    rxn = cobra.Reaction(rxnid)
    model.add_reactions([rxn])
    rxn.reaction = eqn
    # add compartment to metabolites, if not already there
    for met in model.reactions.get_by_id(rxnid).metabolites.keys():
        if type(met.compartment) == type(None):
            met.compartment = met.id.split('_')[-1]
    rxn.bounds = (0,0)

    # add to RBA model
    met_dict = metabolites_dict_from_reaction_equation_RBA(rxn.reaction)
    for k,v in met_dict.items():
        if v != '' and k != '':
            sij_add.add("'MET-" + k + "'.'" + df_rxns.id[i] + "' " + str(v))
            new_species.add("'MET-" + k + "'")

    ## Example of how you might count 2 related products together
    if p == '3hpp':
        prod_rxns['3hppa'][df_rxns.id[i]] = rxnid
        prod_rxns['3hppb'][df_rxns.id[i]] = rxnid
    prods.add(p)
    if p in prod_rxns.keys():
        prod_rxns[p][df_rxns.id[i]] = rxnid
    else:
        prod_rxns[p] = dict()
        prod_rxns[p][df_rxns.id[i]] = rxnid
    added_rxns.add(df_rxns.id[i])
    added_rxns_dict[df_rxns.id[i]] = rxnid
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
            # if the cell is empty or NaN, set it to ''
            if pd.isna(df_prod_info.loc[i,j]):
               prod_info[df_prod_info.loc[i,0]][colnames[j-1]] = ''
            prod_info[df_prod_info.loc[i,0]][colnames[j-1]] = df_prod_info.loc[i,j]
            # can add simulation settings for each product (from SimulationSettings column), to reflect conditions used to make it (where applicable)

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

cobra.io.save_json_model(model, newmodel_path)

# %% [markdown]
#  # Run FBA

if rerun_FBA:
    fba_results = {p:dict() for p in prods}

    for p in prods:
        with model as m:
            while True: # to allow for standardized concluding lines and break points at specific spots
                # write to simulation options JSON file
                fba_settings = dict()
                # fba_results[p] = {'disclaimers':[]} # to preserve order of disclaimers
                # biomassRxn = find_biomass_reactions(m)[0]
                biomassRxn = m.reactions.get_by_id(FBA_biomass_rxn)
                c_sources = []
                if p in ['3hppa', '3hppb','3hpp']:
                    objrxns = 'EX_3hpp_e'
                    p_for_mw = '3hppa'
                else:
                    objrxns = 'EX_' + p + '_e'
                    p_for_mw = p
                # create folder in outputFolder for this product if not there already
                if not os.path.exists(outputFolder + p):
                    # copy from (template) folder
                    shutil.copytree(outputFolder + '(template)', outputFolder + p)
                if p in ['citm', 'polyhb']:
                    m.reactions.THRA_c.knock_out()
                fba_settings['vprod_coreid'] = objrxns
                fba_settings['vprod'] = 'RXN-'+objrxns+'_FWD-SPONT' # default if not found (e.g., due to already being in model)
                for k,v in prod_rxns[p].items():
                    m.reactions.get_by_id(v).bounds = (0,1000)
                    if v == objrxns:
                        fba_settings['vprod'] = k
                # default settings
                # model media - same as scRBA
                # model.reactions.ATPM_c.lower_bound = 1
                c_sources.append({'rxn':'EX_glc__D_e','MW':180.156})
                m.reactions.EX_glc__D_e.bounds = (-13.21,1000)
                biomassRxn.bounds = (0.1, 0.1)
                m.reactions.EX_o2_e.bounds = (-26.25403530509562,1000)
                m.reactions.EX_co2_e.bounds = (0,1000)

                # Allow YNB media uptake
                m.reactions.EX_thm_e.bounds = (-1000,1000)
                m.reactions.EX_ribflv_e.bounds = (-1000,1000)
                m.reactions.EX_nac_e.bounds = (-1000,1000)
                m.reactions.EX_pydxn_e.bounds = (-1000,1000)
                m.reactions.EX_fol_e.bounds = (-1000,1000)
                m.reactions.EX_pnto__R_e.bounds = (-1000,1000)
                m.reactions.EX_4abz_e.bounds = (-1000,1000)
                m.reactions.EX_inost_e.bounds = (-1000,1000)
                # implement simulation settings if they're not empty
                if prod_info[p_for_mw]['SimulationSettings'] not in ['','None','NaN'] and type(prod_info[p_for_mw]['SimulationSettings']) == str:
                    simulation_settings = prod_info[p_for_mw]['SimulationSettings']
                    # contents are initially a string with semicolon-separated key-value pairs; make it a dict
                    simulation_settings = simulation_settings.split(';')
                    simulation_settings = [x.split(':') for x in simulation_settings]
                    for k,v in simulation_settings:
                        if k == 'media':
                            if v == 'reset':
                                # set all uptake bounds to 0
                                for rxn in m.reactions:
                                    if 'EX_' in rxn.id:
                                        rxn.lower_bound = 0
                        elif k == 'c_sources':
                            # set carbon sources
                            c_sources = []
                            for r in v.split(','):
                                # find MW of carbon source from formula
                                # 1. find metabolite ID
                                for met in metabolites_dict_from_reaction_equation_RBA(m.reactions.get_by_id(r).reaction).keys():
                                    if met.endswith('_e'):
                                        # find MW
                                        mass = m.metabolites.get_by_id(met).formula_weight
                                        c_sources.append({'rxn':r,'MW':mass})
                        elif m.reactions.has_id(k):
                            # if a tuple is given, remove the parentheses and split by comma
                            if '(' in v:
                                v = v.replace('(','').replace(')','')
                                v = v.split(',')
                            m.reactions.get_by_id(k).bounds = (float(v[0]),float(v[1]))
            # check if RBA uptakes/biomass differ (i.e., if 'new_FBA_constraints.csv' exists) and apply them
                if os.path.exists(outputFolder + p + '/new_FBA_constraints.csv'):
                    print('Applying RBA uptakes and growth rate to FBA')
                    with open(outputFolder + p + '/new_FBA_constraints.csv', 'r') as f:
                        for line in f:
                            rxn,v = line.split(',')
                            print('old bounds: '+ rxn + ': ' + str(m.reactions.get_by_id(rxn).bounds))
                            # if rxn.startswith('EX_'):
                            #     # m.reactions.get_by_id(rxn).lower_bound = float(v)
                            # else:
                            #     m.reactions.get_by_id(rxn).bounds = (float(v),float(v))
                            print('bound changed: '+ rxn + ': ' + str(m.reactions.get_by_id(rxn).bounds))
                # objective: maximize yields by maximizing objrxns and minimizing substrate uptake
                # objrxns coefficient = MW of product
                m.reactions.get_by_id(objrxns).objective_coefficient = 1
                # m.reactions.get_by_id(objrxns).objective_coefficient = prod_info[p_for_mw]['MW (g/mol)']
                # for met in c_sources:
                #     m.reactions.get_by_id(met['rxn']).objective_coefficient = -met['MW']
                #     print('Objective: ' + met['rxn'] + ' (' + str(m.reactions.get_by_id(met['rxn']).objective_coefficient) + ')')
                # print('Objective: ' + objrxns + ' (' + str(m.reactions.get_by_id(objrxns).objective_coefficient) + ')')
                
                # linear fractional programming for yield maximization: maximize objrxns / -(sum of fluxes for all c_sources)
                # m.objective = m.problem.Objective(m.reactions.get_by_id(objrxns).flux_expression / -sum([m.reactions.get_by_id(x['rxn']).flux_expression for x in c_sources]), direction='max')
                # m.objective = m.problem.Objective(prod_info[p_for_mw]['MW (g/mol)']*m.reactions.get_by_id(objrxns).flux_expression / -sum([x['MW']*m.reactions.get_by_id(x['rxn']).flux_expression for x in c_sources]), direction='max')

                # if exp. yields are provided, set bounds to them
                # add constraint enforcing minimum yields: flux(objrxns)*(product MW) ≥ flux(c_source1)*(c_source1 MW) + flux(c_source2)*(c_source2 MW) + ...
                if pd.notnull(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']):
                    # minYield = m.problem.Constraint(prod_info[p_for_mw]*m.reactions.get_by_id(objrxns).flux_expression - sum([prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']*m.reactions.get_by_id(x['rxn']).flux_expression*x['MW'] for x in c_sources]),lb=0,ub=vmax)
                    minYield = m.problem.Constraint(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']*m.reactions.get_by_id(objrxns).flux_expression - sum([x['MW']*m.reactions.get_by_id(x['rxn']).flux_expression for x in c_sources]), lb=0, ub=vmax)
                    m.add_cons_vars(minYield)
                # if pd.notnull(prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))']) and pd.notnull(prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))']) and prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'] != prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))']:
                #     m.reactions.get_by_id(objrxns).bounds = (prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'],prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))'])
                # biomassRxn.bounds = (0,1000)
                # m.reactions.ATPM_c.bounds = (0,1000)
                # adp = model.metabolites.adp_c
                # gam = biomassRxn.metabolites[adp]
                # # remove gam from biomass
                # for met in ['atp_c','h2o_c']:
                #     biomassRxn.add_metabolites({m.metabolites.get_by_id(met):gam})
                # for met in ['adp_c','h_c','pi_c']:
                #     biomassRxn.add_metabolites({m.metabolites.get_by_id(met):-gam})
                fba = m.optimize()
                # store uptake bounds in case they change later
                uptake_bounds_initial = {x.id:x.lower_bound for x in m.reactions if x.id.startswith('EX_') and x.lower_bound < 0}
                if fba.status != 'optimal':
                    msg = 'FBA for ' + p + ' failed'
                    print(msg)
                    # fba_results[p]['disclaimers'].append(msg)
                    if not adjust_constraints_if_infeas:
                        break
                    biomass_bounds = biomassRxn.bounds
                    uptake_bounds = {x.id:x.lower_bound for x in m.reactions if x.id.startswith('EX_') and x.lower_bound < 0}
                    uptake_mets = {x:x.replace('EX_','') for x in uptake_bounds.keys()}
                    
                    # find minimum total mass of carbon sources to support growth
                    if c_sources == []:
                        print('No carbon sources provided')
                    else:
                        # set objrxns bounds to exp. yields
                        if pd.notnull(prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))']) and pd.notnull(prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))']) and prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'] != prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))']:
                            # find highest growth rate possible at exp. yields and uptake rates
                            print('Finding max growth rate for ' + p)
                            biomassRxn.bounds = (0,1000)
                            m.objective = biomassRxn
                            m.reactions.get_by_id(objrxns).bounds = (prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'],prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))'])
                            fba = m.optimize()
                            if fba.status != 'optimal':
                                print('FBA finding max growth rate for ' + p + ' failed')
                            else:
                                print('Max growth rate: ' + str(fba.objective_value))
                                biomassRxn.bounds = (fba.objective_value,fba.objective_value)
                        else:
                            # find the lowest C uptake needed to support experimental yields under the given conditions
                            # make objective function the sum of all carbon sources * their millimolar masses
                            m.objective = {m.reactions.get_by_id(x['rxn']):x['MW']/1000 for x in c_sources}
                            print(m.objective)
                            # remove bounds on c_sources
                            for x in c_sources:
                                m.reactions.get_by_id(x['rxn']).bounds = (-1000,1000)
                                print(x['rxn'],m.reactions.get_by_id(x['rxn']).bounds)
                            # set bounds of objrxns to exp. yields, to see minimum uptake needed to achieve it
                            if pd.notnull(prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))']) and pd.notnull(prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))']) and prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'] != prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))']:
                                m.reactions.get_by_id(objrxns).bounds = (prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'],vmax)
                            # if len({prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'],prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))']}.intersection({'','None',None,'NaN','nan'})) < 1:
                            #     m.reactions.get_by_id(objrxns).bounds = (prod_info[p_for_mw]['fluxExp_LB_(mmol product/(gDW*h))'],prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))'])
                            else:
                                print('No experimental flux bounds provided for ' + p + '; finding max growth instead')
                                # remove biomass bounds
                                biomassRxn.bounds = (0,1000)
                                # maximize growth
                                m.objective = biomassRxn
                                print(m.objective)
                                fba = m.optimize()
                                if fba.status != 'optimal':
                                    print('FBA finding max growth via ' + objrxns + ' failed')
                                else:
                                    print('Max growth rate: ' + str(fba.objective_value))
                                # print('Uptake bounds: ' + str(uptake_bounds))
                                # print('Uptake metabolites: ' + str(uptake_mets))
                                for i in uptake_bounds:
                                    m.reactions.get_by_id(i).bounds = (-1000,1000)
                                    # minimize difference between flux and old bounds, adding a penalty for each unit of flux
                                # find max growth
                                m.objective = biomassRxn
                                fba = m.optimize()
                                # compare to biomass_bounds (lower)
                                max_growth = fba.objective_value
                                if max_growth < biomass_bounds[0]:
                                    print('When removing bounds on uptakes, max growth rate: ' + str(max_growth) + ' (lower than ' + str(biomass_bounds[0]) + ')')
                                    # force growth to max value
                                    biomassRxn.lower_bound = max_growth
                                else:
                                    print('Max growth rate: ' + str(max_growth))
                                    # force growth to be at least as high as the lower bound
                                    biomassRxn.bounds = (biomass_bounds[0],1000)
                                # find lowest uptakes (mass-adjusted) that support this growth rate
                                m.objective = {m.reactions.get_by_id(k):m.metabolites.get_by_id(v).formula_weight/1000 for k,v in uptake_mets.items()}
                                fba = m.optimize()
                                if fba.status != 'optimal':
                                    print('FBA finding minimum uptake masses for max growth failed')
                                else:
                                    print('Minimum total mass flux (g/gDW/h) needed to support max growth: ' + str(fba.objective_value))
                                    print('Flux for each source:')
                                    # print all non-zero uptake fluxes
                                    for rxn in uptake_mets:
                                        if fba[rxn] != 0:
                                            if uptake_bounds['EX_'+uptake_mets[rxn]] > fba[rxn]:
                                                print("uptake exceeds normal bounds: "+uptake_mets[rxn] + ': ' + str(fba[rxn]))
                                            else:
                                                print(uptake_mets[rxn] + ': ' + str(fba[rxn]))
                                    # print objrxn flux
                                    print(objrxns + ': ' + str(fba[objrxns]))
                                    # store exchange fluxes in CSV file
                                    with open(outputFolder + p + '/fba_exchange_fluxes.csv', 'w') as f:
                                        for rxn in fba.fluxes.index:
                                            if fba[rxn] != 0 and rxn.startswith('EX_'):
                                                f.write(rxn + ',' + str(fba[rxn]) + '\n')
                                        # write growth rate too
                                        f.write('growth,' + str(fba[biomassRxn.id]) + '\n')
                                    break
                            with m as m_test:
                                # remove upper bounds on biomass, maximize it
                                biomassRxn.upper_bound = vmax
                                m_test.objective = biomassRxn
                                fba = m_test.optimize()
                                # use o2 uptake from this simulation as upper bound in the future, to help avoid unrealistic uptakes
                                if fba.status == 'optimal':
                                    m.reactions.EX_o2_e.upper_bound = fba['EX_o2_e']
                            fba = m.optimize()
                            # find lowest uptakes (mass-adjusted) that support this growth rate
                            if fba.status != 'optimal':
                                for i in uptake_bounds:
                                    if p in ['odealc','pmtalc','stalc']: # specific to iRhto models
                                        if i not in ['EX_asp__L_e','EX_thr__L_e','EX_ser__L_e','EX_glu__L_e','EX_gly_e','EX_tyr__L_e','EX_lys__L_e','EX_arg__L_e','EX_trp__L_e','EX_chol_e','EX_ala__L_e','EX_val__L_e','EX_met__L_e','EX_ile__L_e','EX_leu__L_e','EX_phe__L_e','EX_his__L_e','EX_pro__L_e']:
                                            m.reactions.get_by_id(i).bounds = (-1000,1000)
                                    else:
                                        m.reactions.get_by_id(i).bounds = (-1000,1000)
                                    # minimize difference between flux and old bounds, adding a penalty for each unit of flux
                                # find max growth
                                m.objective = biomassRxn
                                fba = m.optimize()
                                # compare to biomass_bounds (lower)
                                max_growth = fba.objective_value
                                if max_growth < biomass_bounds[0]:
                                    print('When removing bounds on uptakes, max growth rate: ' + str(max_growth) + ' (lower than ' + str(biomass_bounds[0]) + ')')
                                    # force growth to max value
                                    biomassRxn.lower_bound = max_growth
                                else:
                                    print('Max growth rate: ' + str(max_growth))
                                    # force growth to be at least as high as the lower bound
                                    biomassRxn.bounds = (biomass_bounds[0],1000)
                                m.objective = {m.reactions.get_by_id(k):m.metabolites.get_by_id(v).formula_weight/1000 for k,v in uptake_mets.items()}
                                # limit upper bounds of uptake fluxes to 0, so solutions with arbitrarily high uptakes offset by high excretion are not found
                                for rxn in uptake_mets:
                                    m.reactions.get_by_id(rxn).upper_bound = 0
                                fba = m.optimize()
                                if fba.status != 'optimal':
                                    print('FBA finding minimum uptake masses for max growth failed')
                                else:
                                    print('Minimum total mass flux (g/gDW/h) needed to support max growth: ' + str(fba.objective_value))
                                    print('Flux for each source:')
                                    # print all non-zero uptake fluxes
                                    for rxn in uptake_mets:
                                        if fba[rxn] != 0:
                                            if uptake_bounds['EX_'+uptake_mets[rxn]] > fba[rxn]:
                                                print("uptake exceeds normal bounds: "+uptake_mets[rxn] + ': ' + str(fba[rxn]))
                                            else:
                                                print(uptake_mets[rxn] + ': ' + str(fba[rxn]))
                                    # print objrxn flux
                                    print(objrxns + ': ' + str(fba[objrxns]))
                                    # store exchange fluxes in CSV file
                                    with open(outputFolder + p + '/fba_exchange_fluxes.csv', 'w') as f:
                                        for rxn in fba.fluxes.index:
                                            if fba[rxn] != 0 and rxn.startswith('EX_'):
                                                f.write(rxn + ',' + str(fba[rxn]) + '\n')
                                        # write growth rate too
                                        f.write('growth,' + str(fba[biomassRxn.id]) + '\n')
                                break
                            if fba.status != 'optimal':
                                print('FBA finding minimum C source uptake masses for ' + objrxns + ' failed')
                                # try again but with no production needed
                                m.reactions.get_by_id(objrxns).bounds = (0,prod_info[p_for_mw]['fluxExp_UB_(mmol product/(gDW*h))'])
                                fba = m.optimize()
                                if fba.status != 'optimal':
                                    print('FBA finding minimum C source uptake masses for ' + objrxns + ' failed when not requiring target product production')
                                    # set glc uptake to 10, maximize biomass
                                    # m.reactions.EX_glc__D_e.bounds = (-10,1000)
                                    # remove biomass bounds
                                    biomassRxn.bounds = (0,1000)
                                    m.objective = biomassRxn
                                    # remove component from biomass rxn and see if it can grow
                                    # # make biomass rxn copy 1st
                                    # m.add_reaction(biomassRxn.copy())
                                    # for met in biomassRxn.reactants:
                                    #     if met.id.startswith('MET-'):
                                    #         # remove component from biomass rxn (- of its coeff)
                                    #         biomassRxn.remove_metabolites(met)
                                    # create demand rxns for each biomass reactant
                                    biomassReactants = [met.id for met in biomassRxn.reactants]
                                    for met in biomassRxn.reactants:
                                        rxn = cobra.Reaction('DM_'+met.id)
                                        m.add_reactions([rxn])
                                        rxn.add_metabolites({met:-1})
                                        # Limit to lower bound of biomass so model is encouraged to produce all parts
                                        rxn.bounds = (0,biomass_bounds[0])
                                    # add to objective
                                    m.objective = {m.reactions.get_by_id('DM_'+k):1 for k in biomassReactants}
                                    fba = m.optimize()
                                    if fba.status != 'optimal':
                                        print('FBA finding minimum C source uptake masses for ' + objrxns + ' failed when not requiring target product production or max biomass')
                                        # try again in YPD media

                                        # Allow YNB media uptake
                                        m.reactions.EX_thm_e.bounds = (-1000,1000)
                                        m.reactions.EX_ribflv_e.bounds = (-1000,1000)
                                        m.reactions.EX_nac_e.bounds = (-1000,1000)
                                        m.reactions.EX_pydxn_e.bounds = (-1000,1000)
                                        m.reactions.EX_fol_e.bounds = (-1000,1000)
                                        m.reactions.EX_pnto__R_e.bounds = (-1000,1000)
                                        m.reactions.EX_4abz_e.bounds = (-1000,1000)
                                        m.reactions.EX_inost_e.bounds = (-1000,1000)
                                        # maximize growth
                                        m.objective = biomassRxn
                                        fba = m.optimize()
                                        if fba.status != 'optimal':
                                            print('FBA finding minimum C source uptake masses for ' + objrxns + ' failed when not requiring target product production or max biomass in YPD media')
                                            # show medium
                                            for rxn in m.reactions:
                                                if rxn.id.startswith('EX_') and rxn.lower_bound < 0:
                                                    print(rxn.id + ': ' + str(rxn.bounds))
                                        else:
                                            print('Max biomass: ' + str(fba.objective_value))
                                    else:
                                        # display missing components
                                        print('Infeasible components for growth: ' + str([k for k in biomassReactants if fba['DM_'+k] == 0]))
                            else:
                                # print objective value and fluxes of all C sources + their names
                                print('Minimum total mass flux (g/gDW/h) needed to match experimental results with target of ' + objrxns + ': ' + str(fba.objective_value))
                                print('Flux for each source:')
                                for rxn in uptake_mets:
                                    if fba[rxn] != 0:
                                        if uptake_bounds['EX_'+uptake_mets[rxn]] > fba[rxn]:
                                            print("uptake exceeds normal bounds: "+uptake_mets[rxn] + ': ' + str(fba[rxn]))
                                            # set uptake bounds to new flux
                                            m.reactions.get_by_id('EX_'+uptake_mets[rxn]).bounds = (fba[rxn],1000)
                                        else:
                                            print(uptake_mets[rxn] + ': ' + str(fba[rxn]))
                                ## find max yield w/ biomass off
                                #biomassRxn.bounds = (0,0)
                                ## turn bounds off for objrxns
                                #m.reactions.get_by_id(objrxns).bounds = (0,1000)
                                #m.objective = m.reactions.get_by_id(objrxns)
                                #fba = m.optimize()
                                #if fba.status == 'optimal':
                                #    print('Max flux w/ biomass off: ' + str(fba.objective_value))
                                #else:
                                #    print('FBA for max flux with biomass off failed')
                            # print all active exchange reactions and their fluxes
                            for rxn in fba.fluxes.index:
                                if fba[rxn] != 0 and rxn.startswith('EX_'):
                                    print('    '+rxn + ': ' + str(fba[rxn]) + ' bounds: ' + str(m.reactions.get_by_id(rxn).bounds))
                    break
                # print all active exchange reactions and their fluxes
                # for rxn in fba.fluxes.index:
                #     if fba[rxn] != 0 and rxn.startswith('EX_'):
                #     # if p == 'tal' and rxn.startswith('EX_') and model.reactions.get_by_id(rxn).bounds != (0,0):
                #         print('    '+rxn + ': ' + str(fba[rxn]) + ' bounds: ' + str(m.reactions.get_by_id(rxn).bounds))
                # find potential carbon source uptake (if not provided already by user); should have negative flux and carbon in its formula
                break
            for n in range(2): # repeat once, in case retries are needed with different settings
                if c_sources == []:
                    for rxn in fba.fluxes.index:
                        if fba[rxn] < 0 and rxn.startswith('EX_'):
                            # find all molecules in the reaction
                            mets = metabolites_dict_from_reaction_equation_RBA(m.reactions.get_by_id(rxn).reaction)
                            for k,v in mets.items():
                                if k.endswith('_e'): # only consider extracellular metabolites (e.g., EX_btn_c should be ignored)
                                    if 'C' in compile_elements_from_formula(m.metabolites.get_by_id(k).formula):
                                        c_sources.append({'rxn':rxn,'MW':m.metabolites.get_by_id(k).formula_weight})
                for c in c_sources:
                    c['v'] = fba[c['rxn']]
                vcarbon = sum([fba[x] for x in [y['rxn'] for y in c_sources]])
                fba_results[p]['vprod'] = fba[objrxns]
                fba_results[p]['vCarbonSources'] = vcarbon
                g_all_substrates = sum([x['v']*x['MW'] for x in c_sources])
                # if more than one carbon source, calculate yield for each
                fba_results[p]['yield'] = abs(prod_info[p_for_mw]['MW (g/mol)']*fba[objrxns]/(g_all_substrates))
                if pd.notnull(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']):
                    sim_to_exp_yield_ratio = fba_results[p]['yield']/prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']
                    break
                    # if less than experimental yields, then run again with altered uptakes
                    # if sim_to_exp_yield_ratio < 1:
                        # # make objective function the sum of all carbon sources * their millimolar masses
                        # m.objective = {m.reactions.get_by_id(x['rxn']):x['MW']/1000 for x in c_sources}
                        # remove bounds on c_sources
                        # for x in c_sources:
                        #     m.reactions.get_by_id(x['rxn']).bounds = (-vmax,vmax)
                        # # set bounds of objrxns based on exp. yields, to see minimum uptake needed to achieve it
                        # m.reactions.get_by_id(objrxns).bounds = (fba[objrxns]*prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']/fba_results[p]['yield'],vmax)

                else:
                    break
            print(objrxns+' flux: '+str(fba[objrxns]))
            for c in c_sources:
                # show yields per g C source, in case it helps
                print(objrxns+' flux / '+c['rxn']+' flux: '+str(fba[objrxns]/c['v']))
            # compare yields to experimental yields (if available)
            print("Simulated yields for " + p + ": " +str(fba_results[p]['yield']) + " g/g all C sources")
            # output GAMS settings (uptake bounds, growth rate, etc.) to file
            if True:
                # add exchange rxn name, MW as keys in list of dicts
                fba_settings['c_sources_list'] = c_sources
                fba_settings['prod_mw'] = prod_info[p_for_mw]['MW (g/mol)']
                # write FBA settings to .py file for use in RBA
                with open(outputFolder + p + '/RBA_defaults_from_FBA.py', 'w') as f:
                    f.write('path_gams = \'../../../GAMS/\'\n')
                    f.write('path_pycore = \'../../../pycore/\'\n')
                    f.write('path_enz_mw = \'../../../input/enz_mw_g_per_mmol.txt\'\n')
                    f.write('path_pro_mw = \'../../../input/pro_mw_g_per_mmol.txt\'\n')
                    f.write('biom_id = \'' + RBA_biomass_rxn + '\'\n')
                    f.write('biom_id_fba = \'' + biomassRxn.id + '\'\n')
                    # write each key from fba_settings as a variable
                    for k,v in fba_settings.items():
                        if type(v) == str:
                            f.write(k + ' = \'' + v + '\'\n')
                        else:
                            f.write(k + ' = ' + str(v) + '\n')
                # write fba_settings to file
                initial_fba_settings = []
                with (open(outputFolder + p + '/RBA_GAMS_defaults_from_FBA.txt', 'w') as f,
                open(outputFolder + p + '/RBA_GAMS_defaults_from_FBA_initial.txt', 'w') as f2, 
                open(outputFolder + p + '/yield_denominators.txt', 'w') as yieldfile,
                open(outputFolder + p + '/yield_denominators_expanded.txt', 'w') as expanded_yieldfile):
                    yieldfile_set = set()
                    expanded_yieldfile_set = set()
                    # list of other possible C sources that may be used if media needs to be updated
                    alternative_c_sources = ['EX_btn_c','EX_glc__D_e','EX_fru_e','EX_glyc_e','EX_lac__D_e','EX_malt_e','EX_sucr_e','EX_xyl__D_e']
                    f.write('$include ../../input/phenotype_common.txt\n') # universal defaults
                    # set vprod to product of interest
                    f.write("$setGlobal vprod '" + fba_settings['vprod'] + "'\n")
                    f.write("parameter minProductYield /"+str(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)'] if pd.notnull(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']) else 0)+"/;\n")
                    f.write("$setGlobal prod_mw "+str(prod_info[p_for_mw]['MW (g/mol)'])+"\n")
                    for rxn in m.reactions:
                        add_to_RBA_bounds = False
                        if rxn.id in added_rxns_dict.values():
                            add_to_RBA_bounds = True
                            for k,v in added_rxns_dict.items():
                                if v == rxn.id:
                                    rba_names = {'fwd':k,'rev':k}
                        elif rxn.id.startswith('RXNADD-'):
                            add_to_RBA_bounds = True
                            rba_names = {'fwd':rxn.id,'rev':rxn.id}
                        elif rxn.id.startswith('EX_'):
                            add_to_RBA_bounds = True
                            rba_names = {'fwd':'RXN-' + rxn.id + '_FWD-SPONT','rev':'RXN-' + rxn.id + '_REV-SPONT'}
                            # for rxns in alternative_c_sources but not in c_sources, add to yieldfile_set
                            for r in alternative_c_sources:
                                if r == rxn.id and r not in [c['rxn'] for c in c_sources]:
                                    expanded_yieldfile_set.add("'"+rba_names['rev']+"' "+str(sum([i.formula_weight for i in m.reactions.get_by_id(r).metabolites])))
                        elif rxn.id == biomassRxn.id:
                            add_to_RBA_bounds = True
                            rba_names = {'fwd':RBA_biomass_rxn,'rev':RBA_biomass_rxn}
                            f.write("$setGlobal mu " + str(fba[biomassRxn.id]) + "\n")
                            f.write("$setGlobal biom_id '" + RBA_biomass_rxn + "'\n")
                            f.write("v.fx('%biom_id%') = %mu% * %nscale%;\n")
                            fba_settings['biom_id'] = RBA_biomass_rxn
                            continue
                        # include bounds in file
                        if add_to_RBA_bounds:
                            file_LB_dict = {f:rxn.lower_bound}
                            if rxn.id in uptake_bounds_initial:
                                file_LB_dict[f2] = uptake_bounds_initial[rxn.id]
                            for file,v_LB in file_LB_dict.items():
                                if rxn.lower_bound == rxn.upper_bound:
                                    file.write("v.fx('" + rba_names['fwd'] + "') = " + str(v_LB) + " * %nscale%;\n")
                                else:
                                    if rxn.lower_bound < 0:
                                        vUptakeRBA = str(-v_LB) # default uptake flux; add carbon slack variable if it's a C source, to use only when RBA is infeasible
                                        for c in c_sources:
                                            if c['rxn'] == rxn.id:
                                                vUptakeRBA = "(" + str(-v_LB) + " + %carbonSlack%)"
                                                yieldfile_set.add("'"+rba_names['rev']+"' "+str(c['MW'])) # create set of compounds to use for yield calculations
                                                expanded_yieldfile_set.add("'"+rba_names['rev']+"' "+str(c['MW']))
                                        if rba_names['rev'] != rba_names['fwd']:
                                            file.write("v.up('" + rba_names['rev'] + "') = " + vUptakeRBA + " * %nscale%;\n")
                                            if rxn.upper_bound < 0:
                                                file.write("v.lo('" + rba_names['rev'] + "') = " + str(-rxn.upper_bound) + " * %nscale%;\n")
                                            else:
                                                file.write("v.lo('" + rba_names['rev'] + "') = 0;\n")
                                            # if rxn is reversible, set bounds for fwd too
                                            if rxn.upper_bound > 0:
                                                file.write("v.up('" + rba_names['fwd'] + "') = " + str(rxn.upper_bound) + " * %nscale%;\n")
                                                file.write("v.lo('" + rba_names['fwd'] + "') = 0;\n")
                                    else:
                                        file.write("v.up('" + rba_names['fwd'] + "') = " + str(rxn.upper_bound) + " * %nscale%;\n")
                                        file.write("v.lo('" + rba_names['fwd'] + "') = " + str(v_LB) + " * %nscale%;\n")
                    yieldfile.write('\n'.join(['/'] + list(yieldfile_set) + ['/']))
                    expanded_yieldfile.write('\n'.join(['/'] + list(expanded_yieldfile_set) + ['/']))
            # store fluxes in CSV file
            with open(outputFolder + p + '/fba_fluxes.csv', 'w') as f:
                # add header
                f.write('rxn,flux,lb,ub\n')
                # write rxn name, flux, and lower and upper bounds
                for rxn in fba.fluxes.index:
                    if fba[rxn] != 0:
                        f.write(rxn + ',' + str(fba[rxn]) + ',' + str(m.reactions.get_by_id(rxn).lower_bound) + ',' + str(m.reactions.get_by_id(rxn).upper_bound) + '\n')
            # print simulated and experimental yields
            if not pd.isnull(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']):
                print(str(fba_results[p]['yield']/prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']*100)+'% of LB of experimental yield: '+str(prod_info[p_for_mw]['yieldExp_LB_(g/g_substrate)']))
            if not pd.isnull(prod_info[p_for_mw]['yieldExp_UB_(g/g_substrate)']):
                print(str(fba_results[p]['yield']/prod_info[p_for_mw]['yieldExp_UB_(g/g_substrate)']*100)+'% of UB of experimental yield: '+str(prod_info[p_for_mw]['yieldExp_UB_(g/g_substrate)']))
            # find yield with growth rate at 10% of previous value (to sustain cell population), then off
            if fba_results[p]['yield'] > 0 and find_min_media_yields:
                # remove upper bounds on objrxns
                m.reactions.get_by_id(objrxns).upper_bound = vmax
                # find minimal media (better for industrial use)
                m.objective = biomassRxn
                max_growth = m.slim_optimize()
                # m.medium = cobra.medium.minimal_medium(m,minimize_components=True)
                # for rxn in m.boundary:
                #     try:
                #         print(rxn.id,rxn.compartments)  # Print the first compartment for each reaction
                #     except IndexError:
                #         print(f"Reaction {rxn.id} has invalid compartments")
                # Note: finding minimal media can be time-consuming
                try:
                    m.medium = cobra.medium.minimal_medium(m,minimize_components=True)
                    print('Minimal media: ' + str(m.medium))
                    # export minimal media to file for GAMS
                    with open(outputFolder + p + '/minimal_media.txt', 'w') as f:
                        for k,v in m.medium.items():
                            # find RBA rxn name
                            if k.startswith('EX_'):
                                rba_name = 'RXN-' + k + '_REV-SPONT'
                            else:
                                rba_name = k
                            f.write("v.up('" + rba_name + "')= " + str(v) + ' * %nscale%;\n')
                except:
                    print('Minimal media not found')
                # find growth rate at 10% of exp.
                m.objective = objrxns
                biomassRxn.bounds = (0.1*fba[biomassRxn.id],0.1*fba[biomassRxn.id])
                fba = m.optimize()
                fba_results[p]['yield_10%_growth'] = prod_info[p_for_mw]['MW (g/mol)']*fba[objrxns]/(g_all_substrates)
                fba_results[p]['growth_10%_exp'] = fba[biomassRxn.id]
                if fba.status == 'optimal':
                    print('Yield at 10% of max growth rate: ' + str(fba_results[p]['growth_10%_exp']) + ' g/g all C sources')
                else:
                    print('FBA for 10% of max growth rate failed')
                # find yield with growth rate off
                biomassRxn.bounds = (0,0)
                fba = m.optimize()
                fba_results[p]['yield_no_growth'] = prod_info[p_for_mw]['MW (g/mol)']*fba[objrxns]/(g_all_substrates)
                if fba.status == 'optimal':
                    print('Yield with growth off: ' + str(prod_info[p_for_mw]['MW (g/mol)']*fba[objrxns]/(g_all_substrates)) + ' g/g all C sources')
                else:
                    print('FBA for growth off failed')
            m.objective = dict() # reset objective

if rerun_RBA:
    # check for 'gams' command
    # try 'module load gams'
    os.system('module load gams')
    if shutil.which('gams') is None:
        print('RBA not run; GAMS not found in PATH')
    else:
        print('Running RBA...')
        for prodname in prods:
            if prodname == '3hpp':
                plist = ['3hppa', '3hppb']
            else:
                plist = [prodname]
            for p in plist:
                cmds = ['cd ' + outputFolder + p,
                        'module load gams',
                        'python3 runRBA.py '+str(int(adjust_constraints_if_infeas))]
                os.system('\n'.join(cmds))
    # read RBA results
    RBA_result_dict = dict()
    # make df to store results; index is rxn name
    df_fluxes = pd.DataFrame(columns=['Rxn'])
    for prodname in sorted(prods):
        if prodname == '3hpp':
            plist = ['3hppa', '3hppb']
        else:
            plist = [prodname]
        for p in plist:
            RBA_result_dict[p] = dict()
            # check if report.txt exists
            if not os.path.exists(outputFolder + p + '/report.txt'):
                print('not found: ' + outputFolder + p + '/report.txt')
                # add empty values for all keys
                for k in ['vprod','vCarbonSources','yield']:
                    RBA_result_dict[p][k] = ''
                continue
            # check if flux_escher.csv exists; add fluxes to fluxes dataframe
            flux_path = outputFolder + p + '/flux.escher.csv'
            if os.path.exists(flux_path):
                with open(flux_path) as f:
                    # convert to dataframe (1st row is header)
                    df_flux = pd.read_csv(f)
                    # add fluxes to df_fluxes
                    # make new column in df_fluxes for this file
                    df_fluxes[p] = ''
                    # for each flux in df_flux, add to df_fluxes (row is rxn name, column is product)
                    for i in df_flux.index:
                        # check if rxn is already in df_fluxes
                        if df_flux.loc[i, 'Rxn'] not in df_fluxes['Rxn'].values:
                            # add row to df_fluxes
                            df_fluxes.loc[df_fluxes.shape[0]] = ''
                            # add rxn name to the index
                            df_fluxes.loc[df_fluxes.shape[0] - 1, 'Rxn'] = df_flux.loc[i, 'Rxn']
                        # add flux to df_fluxes
                        df_fluxes.loc[df_fluxes['Rxn'] == df_flux.loc[i, 'Rxn'], p] = df_flux.loc[i, 'Flux']
            with open(outputFolder + p + '/report.txt') as f:
                text = f.read().split('\n')
            # also read runRBA.yields.txt, if it exists
            fba_results[p]['misc'] = dict()
            if os.path.exists(outputFolder + p + '/runRBA.yields.txt'):
                with open(outputFolder + p + '/runRBA.yields.txt') as f:
                    # convert to dataframe (1st row is header)
                    df = pd.read_csv(f, sep='\t')
                    # find max "growth rate" row
                    mu_max = df['growth rate'].max()
                    # if "method" column is "fba", then add to fba_results
                    if 'method' in df.columns and 'yield' in df.columns:
                        for i in df.index:
                            if df.loc[i,'method'] == 'fba':
                                # add key "yield at mu="+str(df.loc[i,'growth rate']) to fba_results[p]['misc']
                                fba_results[p]['misc']['yield on minimal media at '+str(round(df.loc[i,'growth rate']/mu_max,0))+'% exp. growth'] = df.loc[i,'yield']
                            else:
                                # add key "yield at mu="+str(df.loc[i,'growth rate']) to RBA_result_dict[p]
                                RBA_result_dict[p]['yield on minimal media at '+str(round(df.loc[i,'growth rate']/mu_max,0))+'% exp. growth'] = df.loc[i,'yield']

            for t in text[1:]:
                k,v = t.split('\t')
                try:
                    RBA_result_dict[p][k] = float(v)
                except:
                    RBA_result_dict[p][k] = v
                
            print(RBA_result_dict[p]['vprod'])
    df_fluxes_sorted = df_fluxes.sort_values(by='Rxn')
    df_fluxes_sorted.to_csv(outputFolder + 'fluxes.csv', index=False)

# copy compiled_results_template.xlsx to compiled_results.xlsx
# don't put file extension; copy sheet into outputFolder
output_sheet_path = outputFolder + 'compiled_results'
shutil.copy(os.path.join('./compiled_results_template.xlsx'),
            os.path.join(output_sheet_path+'.xlsx'))
# write FBA yields to rateFBA column (2nd how has headers)
# read headers only
df = pd.read_excel(output_sheet_path+'.xlsx', header=1, nrows=0)
c = 0
for prodname in sorted(prods):
    if prodname == '3hpp':
        plist = ['3hppa', '3hppb']
    else:
        plist = [prodname]
    for p in plist:
        c += 1
        # # add each column from prod_info to df
        df.loc[c, 'product'] = p
        for k,v in prod_info[p].items():
            df.loc[c, k] = v
        if 'misc' in fba_results[p].keys():
            for k,v in fba_results[p]['misc'].items():
                # add column if not there
                if k not in df.columns:
                    df[k] = ''
                df.loc[c, 'FBA: '+k] = v
        for k,v in RBA_result_dict[p].items():
            if k not in df.columns:
                df[k] = ''
            df.loc[c, 'RBA: '+k] = v
        # for k,v in {'FBA':fba_results,'RBA':RBA_result_dict}.items():
        #     for col in ['vprod','yield','vCarbonSources']:
        #         df.loc[c, col+k] = v[p][col]
        # add row for each product
        df.loc[c, 'rateFBA'] = fba_results[p]['vprod']
        df.loc[c, 'vCarbonSourcesFBA'] = fba_results[p]['vCarbonSources']
        df.loc[c, 'yieldFBA'] = fba_results[p]['yield']
        df.loc[c, 'rateRBA'] = RBA_result_dict[p]['vprod']
        df.loc[c, 'vCarbonSourcesRBA'] = RBA_result_dict[p]['vCarbonSources']
        df.loc[c, 'yieldRBA'] = RBA_result_dict[p]['yield']
        try:
            df.loc[c, 'yRBA/yFBA'] = RBA_result_dict[p]['yield']/fba_results[p]['yield']
        except:
            continue
df.to_excel(output_sheet_path+'.xlsx', index=False)
df.to_csv(output_sheet_path+'.csv', index=False)
# store fba_results and RBA_result_dict in json files
with open(fba_results_path, 'w') as f:
    json.dump(fba_results, f)
with open(rba_results_path, 'w') as f:
    json.dump(RBA_result_dict, f)
