from RBA_defaults_from_FBA import *

import sys
sys.path.append(path_pycore)

adjust_constraints_if_infeas = False # default
if len(sys.argv) > 1:
    # There are additional arguments
    adjust_constraints_if_infeas = bool(int(sys.argv[1]))

import json
from simulate import get_GAMS_modelStat, RBA_result

import shutil
shutil.copy(path_gams + 'application/runRBA_max_prod.gms', './runRBA_max_prod.gms');
shutil.copy(path_gams + 'application/soplex.opt', './soplex.opt');

import os
# remove report file if it exists, to avoid accidentally reporting old results
report_path = 'report.txt'
if os.path.exists(report_path):
    os.remove(report_path)
import pandas as pd

# Set growth and glucose uptake rates
carbonSlack = 0 # in case the model needs a rerun, helps test if it can grow

# can hide output (e.g., .lst files) by redirecting to o=/dev/null
hide_output = True
output_redirect_str = ' o=/dev/null' if hide_output else ''

# Initiate report
report = {k:None for k in ['stat', 'vCarbonSources', 'carbonSlack', 'vprod', 'yield']}

# Execute GAMS
os.system('module load gams\n' + 'gams runRBA_max_prod.gms' + \
          ' --carbonSlack=' + str(carbonSlack) + \
          ' --adjust_constraints_if_infeas=' + str(int(adjust_constraints_if_infeas)) + \
          output_redirect_str)
stat = get_GAMS_modelStat('./runRBA_max_prod.modelStat.txt')

def total_carbon_mass_flux(RBA_results, carbon_source_list):
    # in g*mmol/(mol*gDW*h)
    # for each c source, multiply its flux by its MW
    # c_source_list: list of dicts with keys: "rxn" and "MW"
    total_mass_flux = 0
    for c_source in carbon_source_list:
        if c_source['rxn'] in RBA_results.metabolic_flux.keys():
            total_mass_flux += RBA_results.metabolic_flux[c_source['rxn']] * c_source['MW']
    return total_mass_flux

optimal = False
if stat == 'infeasible':
    report['stat'] = stat
    report['vprod'] = 0
    report['yield'] = 0
        
elif stat == 'optimal':
    optimal = True
    res = RBA_result(biom_id=biom_id_fba)
    # convert all RXNADD rxns in runRBA_max_prod.flux.txt to RXN-XXX
    flux_file = './runRBA_max_prod.flux.txt'
    new_flux_file = './runRBA_max_prod.fluxes_rxnadd_as_rxn.txt'
    with open(flux_file) as f:
        text = f.read().split('\n')
    text = [i for i in text if i != '']
    for i in range(len(text)):
        if text[i].split('\t')[0].split('-')[0] == 'RXNADD':
            text[i] = text[i].replace('RXNADD-', 'RXN-')
        # replace biom_id with biom_id_fba
        elif text[i].split('\t')[0] == biom_id:
            text[i] = text[i].replace(biom_id, biom_id_fba,1)
    with open(new_flux_file, 'w') as f:
        f.write('\n'.join(text))
    res.load_raw_flux(filepath=new_flux_file)
    res.calculate_metabolic_flux()
    res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod.replace('RXNADD-', 'RXN-')]
    pflux = res.metabolic_flux[vprod_coreid]
        
elif stat == 'need_rerun':
    itermax = 100; iternum = 0;
    while stat == 'need_rerun' and iternum < itermax:
        iternum += 1
        carbonSlack += 1e-3 # test if it grows with extra substrate
        os.system('module load gams\n' + 'gams runRBA_max_prod.gms' + \
                  ' --carbonSlack=' + str(carbonSlack) + \
                  output_redirect_str)
        stat = get_GAMS_modelStat('./runRBA_max_prod.modelStat.txt')
            
        if stat == 'infeasible':
            report['stat'] = stat
            report['vprod'] = 0
            report['yield'] = 0
            
        elif stat == 'optimal':
            optimal = True
            res = RBA_result(biom_id=biom_id)
            res.load_raw_flux(filepath='./runRBA_max_prod.flux.txt')
            res.calculate_metabolic_flux()
            if vprod.split('-')[0] == 'RXNADD':
                res.metabolic_flux[vprod_coreid] = res.raw_flux[vprod]
            pflux = res.metabolic_flux[vprod_coreid]
        elif stat == 'need_rerun':
            report['stat'] = stat
        
        else:
            print('wtf')

report['carbonSlack'] = carbonSlack
if optimal:
    report['stat'] = stat
    report['vprod'] = pflux
    report['yield'] = pflux * prod_mw / abs(total_carbon_mass_flux(res, c_sources_list))
    report['vCarbonSources'] = sum([res.metabolic_flux[c['rxn']] for c in c_sources_list])

    # Write report text file
    text = []
    for k in report.keys():
        text.append(k + '\t' + str(report[k]))
    with open(report_path, 'w') as f:
        f.write('\n'.join(text))

    # Write JSON results
    with open(path_enz_mw) as f:
        text = f.read().split('\n')
    text = [i for i in text if i != '']
    enz_mw = dict()
    for i in text:
        k,v = i.split('\t')
        enz_mw[k] = float(v)
        
    with open(path_pro_mw) as f:
        text = f.read().split('\n')
    text = [i for i in text if i != '']
    pro_mw = dict()
    for i in text:
        k,v = i.split('\t')
        pro_mw[k] = float(v)
        
    res.enzyme_mw = enz_mw
    res.protein_mw = pro_mw

    res.calculate_all()
    res.enzyme_mw = ''
    res.protein_mw = ''

    res.save_to_json('./RBA_result.json')
    res.make_escher_csv('./flux.escher.csv')
    # add biomass and BIOSYN rxns to alternate csv file
    # res.make_escher_csv('./flux_with_biomass.escher.csv',RBA_rxns_to_include=[biom_id_fba],RBA_regex_to_include=[r'^(?!'+f'{biom_id}'+r')(?=BIOSYN-).*'])
    res.make_escher_csv('./flux_with_biomass.escher.csv',RBA_rxns_to_include=[biom_id_fba])
    # find FBA fluxes from fba_fluxes.csv
    fba_fluxes = pd.read_csv('fba_fluxes.csv')
    
    fba_and_rba_identical_uptakes_and_biomass = True
    # if new_FBA_constraints.csv exists, remove it (in case it's from a prior run)
    new_FBA_constraints_name = 'new_FBA_constraints.csv'
    if os.path.exists(new_FBA_constraints_name):
        os.remove(new_FBA_constraints_name)
    # compare FBA and RBA fluxes, to rerun FBA with RBA constraints if needed
    for rxn in res.metabolic_flux.keys():
        if rxn in fba_fluxes['rxn'].values:
            fba_flux = fba_fluxes[fba_fluxes['rxn'] == rxn]['flux'].values[0]
            rba_flux = res.metabolic_flux[rxn]
            # if rxn is uptake or biomass, compare fluxes
            if (rxn[:2] == 'EX' and fba_flux < 0) or rxn == res.biom_id:
                # replace biomass name with biom_id_fba
                rxn_name = rxn
                if rxn == res.biom_id:
                    rxn_name = biom_id_fba
                if abs(fba_flux - rba_flux) > 1e-3:
                    fba_and_rba_identical_uptakes_and_biomass = False
                    print('FBA and RBA fluxes differ for ' + rxn + ': ' + str(fba_flux) + ' vs ' + str(rba_flux))
                    # write new FBA constraints, for use in A1
                    with open(new_FBA_constraints_name, 'a') as f:
                        f.write(rxn_name + ',' + str(rba_flux) + '\n')
