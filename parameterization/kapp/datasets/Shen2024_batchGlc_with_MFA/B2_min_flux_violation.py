# update model-specific settings in kapp_options.py
from kapp_options import *

path_gams = '../../../../GAMS/parameterization/min_flux_violation/'
path_rxns_list = '../../../../GAMS/model/RBA_rxns.txt'
path_out = './min_flux_violation/'
path_enz_level = './enz_from_proteome/enz_flux_calculation.txt'

run_setting_file_from = './GAMS_setting_files/min_flux_violation_GAMS_settings.txt'
run_setting_file_to = './min_flux_violation/min_flux_violation_GAMS_settings.txt'

df_prot = read_spreadsheet(prot_path)
df_prot.index = df_prot.id.to_list()
min_kapps = dict() # for determining allocation of non-measured proteins w/ used rxns

#### Create directory and copy run settings
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
shutil.copy(run_setting_file_from, run_setting_file_to);

data_val = []; data_idx = [];
for i in df_data.index:
    data_idx.append("'PROSYN-" + df_data.id[i] + "'")
    data_val.append("'PROSYN-" + df_data.id[i] + "' " + str(df_data.loc[i, 'vtrans (mmol/gDW/h)']))

data_val = ['/'] + data_val + ['/']
data_idx = ['/'] + data_idx + ['/']
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data_val))
with open(os.path.join(path_out, 'proteome_data_set.txt'), 'w') as f:
    f.write('\n'.join(data_idx))

#### Determine likely active and inactive reactions informed by proteomics data and calculation in B1_enz_from_proteome.py
# List out model reactions
with open(path_rxns_list) as f:
    idx = f.read().split('\n')[1:-1]
idx = [i[1:-1] for i in idx]
rxns_all = [i for i in idx if i.split('-')[0] == 'RXN']
rxns_enz_all = [i for i in idx if i.split('-')[0] == 'ENZLOAD']
rxns_enz_all = ['RXN-'+i[8:] for i in rxns_enz_all]
enz_rxn_pairs = dict()
# define all enz-rxn pairs, for determining prototype kapps
for rxn in rxns_enz_all:
    tag,rxn_base_id,rxn_dir,enz_id = extract_details_from_rxnid(rxn)
    # find all rxns catalyzed by the same enzyme
    if enz_id not in enz_rxn_pairs.keys():
        enz_rxn_pairs[enz_id] = dict()
    enz_rxn_pairs[enz_id][rxn] = {'tag':tag,'rxn_base_id':rxn_base_id,'rxn_dir':rxn_dir,'v':0,'v_enzload':0,'kapp_proto':0}

rxns_nonenz = [i for i in rxns_all if i not in rxns_enz_all]

# List active enzymatic reactions (rxns with enzymes that can be made, using estimates from B1_enz_from_proteome.py)
with open(path_enz_level) as f:
    enz_fluxes = f.read().split('\n')
enz_fluxes = [i for i in enz_fluxes if i != '']

rxns_enz_active = dict()
for i in enz_fluxes:
    enzid,v = i.split('\t')
    etype = enzid.split('-')[0]
    rxn = 'RXN-'+enzid[8:]
    tag,rxn_base_id,rxn_dir,enz_id = extract_details_from_rxnid(rxn)
    if etype == 'ENZLOAD':
        if float(v) > 0:
            rxns_enz_active[rxn] = {'enz':enz_id,'enzload':enzid,'v_enzload':float(v),'v':0}
            enz_rxn_pairs[enz_id][rxn]['v_enzload'] = float(v)

# List of enzymatic reactions with no proteomics data 
rxns_enz_with_no_prodata = [i for i in rxns_all if i not in list(rxns_enz_active.keys()) + rxns_nonenz]
rxns_enz_with_no_prodata = ["'" + i + "'" for i in rxns_enz_with_no_prodata]
enzload_with_no_prodata = [i.replace('RXN-','ENZLOAD-',1) for i in rxns_enz_with_no_prodata]
rxns_enz_with_no_prodata = ['/'] + rxns_enz_with_no_prodata + ['/']
fname = os.path.join(path_out, 'enzload_enz_with_no_prodata.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + enzload_with_no_prodata + ['/']))
fname = os.path.join(path_out, 'rxns_enz_with_no_prodata.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(rxns_enz_with_no_prodata))

#### Simulation
# shutil.copy(os.path.join(path_gams, 'min_flux_violation.gms'),
#             os.path.join(path_out, 'min_flux_violation.gms'));
# shutil.copy(os.path.join(path_gams, 'soplex.opt'),
#             os.path.join(path_out, 'soplex.opt'));
# shutil.copy(os.path.join(path_gams, 'enz_alloc.gms'),
#             os.path.join(path_out, 'enz_alloc.gms'));

cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_violation.gms --nscale=' + str(nscale) + output_redirect_str]
os.system('\n'.join(cmds))
stop_if_run_failed('./min_flux_violation/min_flux_violation.modelStat.txt')
#### Convert GAMS-scaled flux to actual flux
# All fluxes
rxns_used = []
rxns_enz_with_no_prodata_essential = dict()
enzyme_rxns = []
enz_rxn_dict = {}

fname = os.path.join(path_out, 'min_flux_violation.rxns_essential_with_no_prodata_gamsscaled.txt')
f2name = os.path.join(path_out, 'min_flux_violation.enzload_essential_with_no_prodata_gamsscaled.txt')
with open(fname) as f, open(f2name,'w') as f2:
    newlines = []
    for i in f.read().split('\n'):
        if i not in ['','/']:
            newlines.append(i.replace('RXN-','ENZLOAD-',1))
            j,v = i.split('\t')
            rxns_enz_with_no_prodata_essential[j[1:-1]] = float(v) 
        else:
            newlines.append(i)
    f2.write('\n'.join(newlines))
    # print(rxns_enz_with_no_prodata_essential)

fname = os.path.join(path_out, 'min_flux_violation.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    tag,rxn_base_id,rxn_dir,enz = extract_details_from_rxnid(r)
    val = float(val) / nscale
    fluxes_new.append('\t'.join([r, vtype, str(val)]))
    rxns_used.append("'" + r + "'")
    # if rxn (as determined by characters after last "-") isn't SPONT or UNKNOWN
    if enz not in spont_rxn_suffixes and r.split('-')[0] == 'RXN':
        enzload = r.replace('RXN-', 'ENZLOAD-',1)
        if enz not in enz_rxn_dict.keys():
            enz_rxn_dict[enz] = []
        enz_rxn_dict[enz].append(enzload)
        enz_rxn_pairs[enz][r]['v'] = float(val)
        if r in rxns_enz_active.keys():
            rxns_enz_active[r]['v'] = val
        enzyme_rxns.append("'" + enzload + "'")
# find min kapps to construct constraints that get the essential enz rxns w/ no prodata to have kapps as close as possible to them
for j in rxns_enz_active.keys():
    v = rxns_enz_active[j]['v']
    if v > 0:
        min_kapp = mu * v / (float(rxns_enz_active[j]['v_enzload']))
        enz_rxn_pairs[rxns_enz_active[j]['enz']][j]['kapp_proto'] = min_kapp
        min_kapps[j] = min_kapp
kapp_benchmark = float(np.median(list(min_kapps.values())))
kapp_max = float(max(list(min_kapps.values())))
print('max predicted kapp',str(kapp_max))
print('benchmark kapp',str(kapp_benchmark))

enzload_constraints = [] # determining enzyme allocation for next problem
kapp_init = []
c = 0

# find kapps for all rxns that don't have any
for enz in enz_rxn_pairs.keys():
    for rxn in enz_rxn_pairs[enz].keys():
        if enz_rxn_pairs[enz][rxn]['kapp_proto'] == 0:
            rxn_base_id = enz_rxn_pairs[enz][rxn]['rxn_base_id']
            v = enz_rxn_pairs[enz][rxn]['v']
            v_enzload = enz_rxn_pairs[enz][rxn]['v_enzload']
            kapps = []
            # # assign kapp as the sum of all fluxes catalyzed by the enzyme, divided by the enzyme's synthesis flux
            # for enz2 in enz_rxn_pairs.keys():
            #     for rxn2 in enz_rxn_pairs[enz2].keys():
            #         if enz_rxn_pairs[enz2][rxn2]['rxn_base_id'] == rxn_base_id:
            #             if enz_rxn_pairs[enz2][rxn2]['kapp_proto'] > 0:
            #                 kapps.append(enz_rxn_pairs[enz2][rxn2]['kapp_proto'])
                            
            #                 min_kapps[rxn] = float(np.median([enz_rxn_pairs[enz2][rxn2]['kapp_proto']]))
            #                 break
            # kapps = [enz_rxn_pairs[enz][i]['kapp_proto'] for i in enz_rxn_pairs[enz].keys() if enz_rxn_pairs[enz][i]['kapp_proto'] > 0]
            # if len(kapps) > 0:
            #     enz_rxn_pairs[enz][rxn]['kapp_proto'] = float(np.median(kapps))
            #     min_kapps[rxn] = float(np.median(kapps))
            # else:
                # look for kapps for the same rxn in other enzymes
                # kapps = [enz_rxn_pairs[i][j]['kapp_proto'] for i in enz_rxn_pairs.keys() for j in enz_rxn_pairs[i].keys() if j['rxn_base_id'] == enz_rxn_pairs[enz][rxn]['rxn_base_id'] and enz_rxn_pairs[i][j]['kapp_proto'] > 0]
                # if len(kapps) > 0:
                #     enz_rxn_pairs[enz][rxn]['kapp_proto'] = float(np.median(kapps))
                #     min_kapps[rxn] = float(np.median(kapps))
                # else:
            if enz_rxn_pairs[enz][rxn]['kapp_proto'] == 0:
                # assign benchmark kapp
                enz_rxn_pairs[enz][rxn]['kapp_proto'] = kapp_benchmark
        # constraint to encourage equal distribution of enzload amongst all used rxns
        c += 1
        enzsyn = 'ENZSYN-' + enz
        enzload = 'ENZLOAD-'+rxn.split('-',1)[1]
        enzload_constraints.append("Equation EnzLoadConstraint" + str(c) + "; EnzLoadConstraint" + str(c) + ".. " + "EnzLoadSlackPos('" + rxn + "') - EnzLoadSlackNeg('" + rxn + "') + v('" + rxn + "') =e= v('" + enzsyn + "') / " + str(len(enz_rxn_pairs[enz].keys())) + ";")
        kapp_init += [f"'{rxn}' {str(enz_rxn_pairs[enz][rxn]['kapp_proto'])}"]

# for enz,rxns in enz_rxn_dict.items():
#     enzsyn = 'ENZSYN-' + enz
#     for rxn in rxns:
#         c += 1
#         c1 += 1
#         eload = rxn.replace('RXN-','ENZLOAD-',1)
#         j = rxn.replace('ENZLOAD-','RXN-',1)
#         if j in rxns_enz_with_no_prodata_essential.keys():
#             kapp_init.append("Equation kappEst" + str(c1) + "; kappEst" + str(c1) + ".. kappEstSlackPos('" + j + "') - kappEstSlackNeg('" + j + "') + v('" + eload + "') =e= v('" + j + "') / (%mu% * " + str(kapp_benchmark) + ");") 
#         elif rxns_enz_active[j]['v'] > 0:
#             kapp_init.append("Equation kappEst" + str(c1) + "; kappEst" + str(c1) + ".. kappEstSlackPos('" + j + "') - kappEstSlackNeg('" + j + "') + v('" + eload + "') =e= v('" + j + "') / (%mu% * " + str(min_kapps[j]) + ");") 
#         # constraint to encourage equal distribution of enzload amongst all used rxns
#         enzload_constraints.append("Equation EnzLoadConstraint" + str(c) + "; EnzLoadConstraint" + str(c) + ".. " + "EnzLoadSlackPos('" + rxn + "') - EnzLoadSlackNeg('" + rxn + "') + v('" + rxn + "') =e= v('" + enzsyn + "') / " + str(len(rxns)) + ";")
fname = os.path.join(path_out, 'kapp_init.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + kapp_init + ['/']))
fname = os.path.join(path_out, 'enz_alloc_constraints.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(enzload_constraints))
fname = os.path.join(path_out, 'min_flux_violation.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
fname = os.path.join(path_out, 'min_flux_violation.rxns_used.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + rxns_used + ['/']))
# for determining enzyme allocation
fname = os.path.join(path_out, 'min_flux_violation.enzyme_rxns.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(['/'] + enzyme_rxns + ['/']))
   
# Only flux of reactions not associated with expressed enzymes
fname = os.path.join(path_out, 'min_flux_violation.flux_essential_with_no_prodata_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    val = float(val) / nscale
    fluxes_new.append('\t'.join([r, vtype, str(val)]))
fname = os.path.join(path_out, 'min_flux_violation.flux_essential_with_no_prodata.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
# for determining enzyme allocation
cmds = ['cd ' + path_out,
        'module load gams',
        'gams enz_alloc.gms --nscale='+str(nscale2)]
os.system('\n'.join(cmds))
stop_if_run_failed('./min_flux_violation/enz_alloc.modelStat.txt')
# All fluxes
fname = os.path.join(path_out, 'enz_alloc.flux_gamsscaled.txt')
rxns_used = []
enzyme_rxns = []
enz_rxn_dict = {}
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
aa_comp = dict()
biomass_from_RBA = dict()
for i in fluxes:
    r,vtype,val_raw = i.split('\t')
    val = float(val_raw) / nscale2
    fluxes_new.append('\t'.join([r, vtype, str(val)]))
    rxns_used.append("'" + r + "'")
    tag,rxn_base_id,rxn_dir,enz_id = extract_details_from_rxnid(r)
    # if rxn (as determined by characters after last "-") isn't SPONT or UNKNOWN
    if enz not in spont_rxn_suffixes and tag == 'RXN':
        enzload = r.replace('RXN-', 'ENZLOAD-',1)
        if enz not in enz_rxn_dict.keys():
            enz_rxn_dict[enz] = []
        enz_rxn_dict[enz].append(enzload)
        enzyme_rxns.append("'" + enzload + "'")
    elif tag == 'PROSYN':
        pro = r.split('-')[-1]
        # calculate amino acid composition, for use in the GSM model (if different/needed)
        # find the amino acid sequence of the protein
        aa_seq = df_prot.loc[pro, 'sequence']
        for aa in aa_seq:
            # check if it's a letter
            if aa.isalpha():
                # add val to aa_comp, since for every mol of protein made, 1 mol of each aa in its sequence is made
                if aa in aa_comp.keys():
                    # aa_comp[aa]['total flux'] += val_raw
                    aa_comp[aa] += float(val_raw)
                else:
                    # aa_comp[aa] = {'total flux':val_raw}
                    aa_comp[aa] = float(val_raw)
        # Optional per-protein amino acid composition calculations
        # json.dump(aa_comp, open(os.path.join(path_out, 'enz_alloc_aa_comp_'+enz+'.json'), 'w'), indent=4)

# divide aa_comp values by biomass flux to get new GSM coefficients, in case new ones are needed
for aa in aa_comp.keys():
    aa_comp[aa] = aa_comp[aa]/nscale2
    aa_biom_coeff = aa_comp[aa]/(mu)
    # find matching tRNA name in aa_map, and create draft of biomass rxn from it (in case it differs from GSM)
    biomass_from_RBA[aa_map[aa_map['aa_abbv'] == aa.upper()]['tRNA_in'].to_list()[0].split('MET-')[1]] = -aa_biom_coeff
    biomass_from_RBA[aa_map[aa_map['aa_abbv'] == aa.upper()]['tRNA_out'].to_list()[0].split('MET-')[1]] = aa_biom_coeff
    
# calculate mol % of each aa in aa_comp
# for aa in aa_comp.keys():
#     aa_comp[aa]['mol %'] = aa_comp[aa]['total flux']/sum([aa_comp[aa]['total flux'] for aa in aa_comp.keys()])
# store aa composition, biomass details in a file
json.dump(biomass_from_RBA, open(os.path.join(path_out, 'enz_alloc_biomass_from_RBA_partial.json'), 'w'), indent=-1, sort_keys=True)
json.dump(aa_comp, open(os.path.join(path_out, 'enz_alloc_total_aa_fluxes_for_prot_unscaled.json'), 'w'), indent=-1, sort_keys=True)
# find lowest flux out of all vtrans values in df_data and all fluxes in fluxes_new, as well as all predicted ENZLOAD fluxes based on EnzLoadConstraint values
all_nonzero_fluxes = df_data[df_data['vtrans (mmol/gDW/h)'] > 0]['vtrans (mmol/gDW/h)'].to_list()
# all_nonzero_fluxes = df_data[df_data['vtrans (mmol/gDW/h)'] > 0]['vtrans (mmol/gDW/h)'].to_list() + [float(i.split('\t')[2]) for i in fluxes_new if float(i.split('\t')[2]) > 0] \
# + [float(i.split('\t')[2])/len(enz_rxn_dict[i.split('\t')[2]]) for i in fluxes_new if str(i.split('\t')[0]).startswith('ENZSYN-')]
enzloads_predicted = dict()
for i in fluxes_new:
    r,vtype,val = i.split('\t')
    val = float(val) * nscale2
    prefix,enz = r.split('-')[0],r.split('-')[-1]
    if val > 0:
        all_nonzero_fluxes.append(val)
        if enz not in spont_rxn_suffixes:
            if prefix == 'RXN':
                enzload = r.replace('RXN-', 'ENZLOAD-',1)
                if enzload not in enzloads_predicted.keys():
                    enzloads_predicted[enzload] = 0
            elif prefix == 'ENZSYN':
                if enz in enz_rxn_dict.keys():
                    v_enzload_pred = val/len(enz_rxn_dict[enz])
                    # add predicted ENZLOAD flux, in case the solver treated it as 0
                    if v_enzload_pred > 0:
                        all_nonzero_fluxes.append(v_enzload_pred)
                    else:
                        print('Predicted ENZLOAD flux for each rxn using', enz, 'is 0')
                else:
                    print('No ENZLOAD reactions found for', enz)
            elif prefix == 'ENZLOAD':
                enzloads_predicted[r] = val
# check if any ENZLOAD fluxes were predicted to be 0 despite their associated rxn fluxes having nonzero values
for enzload,v in enzloads_predicted.items():
    if v == 0:
        print('Predicted ENZLOAD flux for', enzload, 'is 0')
# all_nonzero_fluxes = list(df_data['vtrans (mmol/gDW/h)']) + [float(i.split('\t')[2]) for i in fluxes_new] \
    # + [float(i.split('\t')[2])/len(enz_rxn_dict[i.split('\t')[2]]) for i in fluxes_new if str(i.split('\t')[0]).startswith('ENZSYN-')]
v_min = min(all_nonzero_fluxes)
v_max = max(all_nonzero_fluxes)
print('v_min:', v_min, 'v_max:', v_max)
print('Make sure both values fit within the solver tolerance range')

fname = os.path.join(path_out, 'enz_alloc.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
