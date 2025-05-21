# update model-specific settings in kapp_options.py
from kapp_options import *
min_flux_cutoff = 0 
# tolerance value added to assist with rounding kapps; set to 0 if not needed
tol = 0 

# Load enzyme info
df_enz = read_spreadsheet('../../../../build_model/model/ENZYME_stoich_curation.tsv')

# Load path of set4 (enzyme-reaction many-to-many mapping)
set4_path = './kapp_ambiguousLoad_case_resolve_common.txt'
# set4 handles cases of multiple enzymes corresponding to multiple rxns (many-to-many)
# columns: 
#   enzyme (for all enz that correspond to a rxn)
#   rxn: 
#       disableManually for any I don't want to count
#       noRxnFlux_resolveByManualCheck: behaves identically to disableManually, but perhaps serves as a warning
#   rule: 
#       divmaxenz for cases where I want to treat the [enzyme] as equaling its maximum value

res_metab = RBA_result(biom_id=biom_id)
flux_path = './min_flux_violation/enz_alloc.flux_gamsscaled.txt'
# check if any enzymatic rxns witho no proteomics data were used
with open('min_flux_violation/enz_alloc.flux_essential_with_no_prodata_unscaled.txt') as f:
    rxns_essential_NP = f.read().split('\n')
rxns_essential_NP = [i for i in rxns_essential_NP if i != '']
rxns_essential_NP = [i.split('\t')[0] for i in rxns_essential_NP]

# extract lines starting with 'ENZSYN from fluxes
with open('./min_flux_violation/enz_alloc.flux_gamsscaled.txt') as f:
    rxns_used = f.read().split('\n')
# find flux of total protein going into biomass: BIOSYN-PROTTOBIO
# make dict of all PROWASTE rxns (except PROWASTE-TOTALPROTEIN) with name as key and flux as value
prowaste_used = dict()
prosyn_used = dict()
enzsyn_used = []
enz_rxns_without_enzload = dict()
for i in rxns_used:
    i_twocol = i.replace('	v	 ','\t').replace('\tv\t','\t')
    if len(i_twocol.split('\t')) < 2:
        print(i_twocol)
        continue
    rxn,v = i_twocol.split('\t')[:2]
    v = float(v)/nscale2
    tag,rxn_base_id,direc,enz = extract_details_from_rxnid(rxn)
    i_twocol = rxn+'\t'+str(v)
    if tag == 'PROWASTE' and rxn != 'PROWASTE-TOTALPROTEIN':
        prowaste_used[rxn] = v 
    elif tag == 'PROSYN':
        prosyn_used[rxn] = v 
    elif tag == 'ENZSYN':
        enzsyn_used.append(i_twocol)
    elif rxn == 'BIOSYN-PROTTOBIO':
        v_prot_to_bio = float(i_twocol.split('\t')[1])
    elif tag == 'RXN' and enz not in spont_rxn_suffixes:
        enz_rxns_without_enzload[rxn] = v

# Add new set of constraints to reflect unused but observed protein production, similar to dummy protein (for every prot, v('PROSYN-prot') >= prowaste_used('PROWASTE-prot') / v_prot_to_bio * %nscale%)

## UNUSED: make separate file for defaults, so not every options file needs to be updated when the model changes
## pro_measured_unused_default_path = '../../../../prosyn_constraints_measured_but_unused_default.txt'
## if os.path.exists(pro_measured_unused_default_path):

pro_measured_unused_path = 'prowaste_from_pro_with_data.txt'
pro_predicted_unused_path = 'prowaste_from_pro_without_data.txt'
c = 0

with open(pro_measured_unused_path,'w') as f, open(pro_predicted_unused_path,'w') as f2, open('prot_mass_flux.txt','w') as f3:
    f3.write('/'+str(v_prot_to_bio)+'/')
    # for each protein in PROWASTE-used, allow flux through PROWASTE-PROT to be as high as infinity (to avoid scaling issues)
    # Also, make the proteome mass fraction comprised of each protein at least equal to the mass fraction comprised of unused copies of it 
    ## (i.e., PROSYN * MW >= ptot * PROWASTE_0 * MW / ptot_0)
    ### divide both sides by MW: PROSYN >= ptot * PROWASTE_0 / ptot_0)
    for prot in prowaste_used.keys():
        prot_id = prot.replace('PROWASTE-','')
        prosyn_id = 'PROSYN-'+prot_id
        if prosyn_id in prosyn_used.keys():
            constraint = "'"+prosyn_id+"' "+str(prowaste_used[prot])+"\n"
            # if protein is in the proteome data (path_data), add to f2
            if prot_id in df_data_full.index:
                f.write(constraint)
            else:
                f2.write(constraint)
            c += 1
# write to new file
kapp_input_enzsyn_path = './kapp_input_enzsyn.txt'
with open(kapp_input_enzsyn_path, 'w') as f:
    f.write('\n'.join(enzsyn_used))
res_metab.load_raw_flux(flux_path,nscale=nscale2)
res_metab.calculate_metabolic_flux()
res_metab.make_escher_csv('./flux_metab.escher.csv')

res_esyn = RBA_result(biom_id=biom_id, twocol_format=True)
# load_raw_flux using a version of min_flux_violation.flux.txt that only contains lines starting with 'ENZSYN
#res_esyn.load_raw_flux('./min_flux_violation/enz_alloc.flux.txt')
# res_esyn.load_raw_flux('./min_flux_violation/enz_flux_calculation.txt')
res_esyn.load_raw_flux(kapp_input_enzsyn_path)
with open(metab_rxns_path,"r") as f: 
    all_rxns = [line.strip().replace("'","") for line in f.readlines() if line.strip() != "/"]
with open(path_model+"SM_rxn_ids.txt","r") as f:  
    all_GSM_rxns = [line.strip().replace("'","") for line in f.readlines() if line.strip() != "/"]

mu = res_metab.growth_rate
print('Growth rate:', mu)
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
    if rxn in rxndict.keys() and enz in spont_rxn_suffixes:
        rxndict[rxn].append('zeroCost')
        
    if rxn in rxndict.keys() and enz in enzdict.keys():
        rxndict[rxn].append(enz)
        enzdict[enz].append(rxn)
        
rxndict = {k:set(v) for k,v in rxndict.items()}
#print('rxndict',rxndict)
rxndict_zeroCost = {k:v for k,v in rxndict.items() if v == {'zeroCost'}}
rxndict = {k:v for k,v in rxndict.items() if 'zeroCost' not in [k,v]}
enzdict = {k:set(v) for k,v in enzdict.items() if 'zeroCost' not in [k,v]}
#print('enzdict',enzdict)

# Find enzymes that together carry a total load of reactions
x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
enz_share_rxn_load = set().union(*[v for v in x.values()])
enz_share_rxn_load.discard('zeroCost')

# Find enzymes that individually carry loads of multiple reactions
x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
enz_multiload = set(x.keys())

# Set 1: Enzyme-reaction one-to-one load mapping
set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])
# write set1 to file
with open('./kapp_set1_one_to_one.txt', 'w') as f:
    f.write('\n'.join(sorted(set1)))

# Set 2: Enzyme-reaction one-to-many
set2 = enz_multiload - enz_share_rxn_load
# write set2 to file
with open('./kapp_set2_one_to_many.txt', 'w') as f:
    f.write('\n'.join(sorted(set2)))

# Set 3: Enzyme-reaction many-to-one
set3 = enz_share_rxn_load - enz_multiload
# write set3 to file
with open('./kapp_set3_many_to_one.txt', 'w') as f:
    f.write('\n'.join(sorted(set3)))

# Set 4: Enzyme-reaction many-to-many
set4 = enz_share_rxn_load & enz_multiload
# write set4 to file
with open('./kapp_set4_many_to_many.txt', 'w') as f:
    f.write('\n'.join(sorted(set4)))
# if set4_path doesn't exist or is only 1 line long, create it with the header
# if not os.path.exists(set4_path) or len(open(set4_path).readlines()) < 2:
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
set4_to_print = set4 - set(enzs) - set(rxns_to_ignore_for_kapps)
if len(set4_to_print) > 0:
    print('set4 enzymes not in set4 file and not ignored:')
    for i in set4_to_print:
        print(i)
#set4 = set(enzs)
    
# From manual checking of the sets, some enzymes in set2 and set3 require special treatment in calculation
# which are recorded in set4. Thus excluded those in set2 and set3
set2 = set2 - set4
set3 = set3 - set4

kapp = dict()
kapps_rba = dict()
### For testing only: set each rxn-enzyme's kapp individually
kapp_minimal_assumptions = dict()
kapp_ma_cutoff = 1e-1
kapp_ma_text = []
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
    try:    
        kapp[rid] = (mu * abs(rval) / enzval) / 3600
        kapps_rba[rid] = kapp[rid] * 3600 + tol
    except:
        print('kapp calc error: mu',type(mu),mu,'rval',type(rval),rval,'enzval',type(enzval),enzval)
    
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
        kapp[rid] = (mu * rvalsum / enzval) / 3600
        kapps_rba[rid] = kapp[rid] * 3600 + tol
        
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
        if enz != 'zeroCost':
            rids.append('RXN-' + rxn + '_' + rdir + '-' + enz)
            enzvals.append(res_esyn.raw_flux['ENZSYN-' + enz])
    enzval = sum(enzvals) # changed to match scRBA suppMat description; originally found the max value, not the sum
        
    for rid in rids:
        kapp[rid] = (mu * abs(rval) / enzval) / 3600
        kapps_rba[rid] = kapp[rid] * 3600 + tol
        
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
            print('no ENZSYN flux found for',rxntext,rule)
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
        #rxn_possible_names = [rxn]
        if rxn not in all_rxns + all_GSM_rxns:
            print('Reaction from '+set4_path+' ignored due to its absence from model: '+rxn)
            continue
        if rxn in res_metab.metabolic_flux.keys():
            rxns_on.append(rxn)
            rxnval += abs(res_metab.metabolic_flux[rxn])
    if rxnval < min_flux_cutoff:
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
        kapp[rid] = (mu * rxnval / enzval) / 3600
        kapps_rba[rid] = kapp[rid] * 3600 + tol
#Flux is numerically low, near zero
path='../exclude_parameterization_list.txt'
if os.path.exists(path):
    with open(path) as f:
        excl = f.read().split('\n')
else:
    excl = []
# # TEST: add all rxns with fluxes below min_flux_cutoff to excl
# for k,v in res_metab.metabolic_flux.items():
#    if abs(v) < min_flux_cutoff:
#        excl.append(k)
excl = [r for r in excl if r != '']
excl += []

texts = []
gsm_kapps = dict()
gsm_kapps_per_hr = dict()
gsm_kapps_with_prodata = dict()
perhr_texts_used_only = []
for k,v in kapp.items():
    if k not in excl:
        texts.append(k + '\t' + str(v))
        gsm_rxn = str(k.replace('RXN-','',1).rsplit('-',1)[-2])
        gsm_rxn_no_dir = gsm_rxn.rsplit('_',1)[-2]
        dir = k.replace('RXN-','',1).rsplit('_',1)[-1]
        if gsm_rxn_no_dir not in all_GSM_rxns:
            print('Reaction from kapp not in GSM: ' + gsm_rxn_no_dir)
        gsm_rxns_listed = list(gsm_kapps.keys())
        if gsm_rxn_no_dir not in gsm_rxns_listed:
            gsm_kapps[gsm_rxn_no_dir] = list()
        gsm_kapps[gsm_rxn_no_dir].append(v)
        if k not in rxns_essential_NP:
            if gsm_rxn_no_dir not in gsm_rxns_listed or gsm_rxn_no_dir not in gsm_kapps_with_prodata.keys():
                gsm_kapps_with_prodata[gsm_rxn_no_dir] = list()
            gsm_kapps_with_prodata[gsm_rxn_no_dir].append(v)

# make RBA versions of gsm_kapps (in 1/h)
for k,v in gsm_kapps.items():
    gsm_kapps_per_hr[k] = [i * 3600 for i in v]

with open('./kapps_in_vivo.txt', 'w') as f:
    f.write('\n'.join(['rxnid\tkapp (1/s)'] + sorted(texts)))
with open('./gsm_kapps.csv', 'w') as f:
    f.write('rxnid,median\n')
    # f.write('rxnid,min,max,median,all_kapps\n')
    # f.write('rxnid,dir,min,max,median,all_kapps\n')
    for k,v in sorted(gsm_kapps.items()):
        cols = [k,str(np.median(v))]
        cols = [k,str(min(v)),str(max(v)),str(np.median(v)),str(v)]
        # cols = [k.rsplit('_',1)[-2],str(k.rsplit('_',1)[1]),str(min(v)),str(max(v)),str(np.median(v)),str(v)]
        f.write(','.join(cols) + '\n')
with open('./gsm_kapps_with_prodata.csv', 'w') as f:
    f.write('rxnid,median\n')
    # f.write('rxnid,min,max,median,all_kapps\n')
    # f.write('rxnid,dir,min,max,median,all_kapps\n')
    for k,v in sorted(gsm_kapps_with_prodata.items()):
        cols = [k,str(np.median(v))]
        cols = [k,str(min(v)),str(max(v)),str(np.median(v)),str(v)]
        # cols = [k.rsplit('_',1)[-2],str(k.rsplit('_',1)[1]),str(min(v)),str(max(v)),str(np.median(v)),str(v)]
        f.write(','.join(cols) + '\n')
kapp_med = np.median(list(kapp.values())) * 3600
kapp_max = max(list(v for v in kapps_rba.values() if v>0)) 
kapp_min = min(list(v for v in kapps_rba.values() if v>0)) 
kapp_ma_med = np.median(list(kapp_minimal_assumptions.values()))
#kapp_ma_max = np.max(list(kapp_minimal_assumptions.values()))
kapp_ma_default = kapp_ma_med
default_kapp = kapp_med + tol
# return kapps in a format suitable for use in RBA
kapp_txt = []

for r in rxns_essential_NP:
    if r not in kapps_rba.keys():
        kapps_rba[r] = default_kapp + tol
        print('No kapp found for ' + r + ', setting to default')
for k,v in kapps_rba.items():
    perhr_texts_used_only.append("'" + k + "'" + '\t' + str(v))

with open('./kapps_per_hr_without_unused_rxns_and_default_kapps.txt', 'w') as f:
    f.write('\n'.join(['/'] + sorted(perhr_texts_used_only) + ['/']))

for k in enz_rxns_without_enzload.keys():
    if k not in kapps_rba.keys():
        kapps_rba[k] = default_kapp + tol
        print('No kapp found for ' + k + ', setting to default')
        perhr_texts_used_only.append("'" + k + "'" + '\t' + str(kapps_rba[k]))

with open('./kapps_per_hr_without_unused_rxns.txt', 'w') as f:
    f.write('\n'.join(['/'] + sorted(perhr_texts_used_only) + ['/']))

output_info = []

# add all other rxns with default kapp
for rxn in all_rxns: 
    tag,rxn_base_id,direc,enz = extract_details_from_rxnid(rxn)
    rxn_name_without_enz = rxn.split('-')[-2]
    if enz not in spont_rxn_suffixes:
        if rxn in kapps_rba.keys():
            if kapps_rba[rxn] <= min_flux_cutoff:
                output_info.append('kapp <= cutoff for ' + rxn)
                kapps_rba[rxn] = default_kapp
        elif rxn in rxns_essential_NP:
            kapps_rba[rxn] = default_kapp # to minimize the risk of kapps making growth infeasible 
        # check if rxn carries flux
        elif rxn_base_id in res_metab.metabolic_flux.keys():
            kapp_found = False
            # check if any isozymes have kapp values in gsm_kapps
            for r in gsm_kapps_per_hr.keys():
                if r == rxn_base_id:
                    kapp_found = True
                    # set to the minimum kapp value for that rxn or the default - whichever is lower.
                    # this reduces the risk of the default kapp being too high and creating unrealistically efficient enzymes
                    output_info.append('Setting kapp of ' + rxn + ' to min kapp of all isozymes')
                    kapps_rba[rxn] = min(min(gsm_kapps_per_hr[r]),default_kapp)
                    # kapps_rba[rxn] = min(gsm_kapps_per_hr[r])
                    break
            if not kapp_found:
                output_info.append('No kapp found for ' + rxn + ' or isozymes of it; setting to default')
                kapps_rba[rxn] = default_kapp
        else: 
            output_info.append('No kapp found for ' + rxn)
            kapps_rba[rxn] = default_kapp
        if rxn not in kapp_minimal_assumptions.keys():
            # set kapp to the max value of all kapps for that enzyme (or reaction, if unavailable), to minimize risk of overestimation of enzyme demand
            new_kapp = 0
            for k,v in kapp_minimal_assumptions.items():
                if enz in k.split('-')[-1]:
                    if v > new_kapp:
                        new_kapp = v + tol
                # find max kapp for that reaction
                elif rxn_name_without_enz == k.split('-')[-2]:
                    if v > new_kapp:
                        new_kapp = v + tol
                # remove location to see if enzyme is identical to another
                elif enz == k.split('-')[-2].split('_')[-1]:
                    if v > new_kapp:
                        new_kapp = v + tol
            if new_kapp == 0:
                output_info.append('No kapp found for ' + rxn + ' or enzyme ' + enz)
                new_kapp = kapp_ma_default
            kapp_minimal_assumptions[rxn] = new_kapp
            kapp_ma_text.append("'" + rxn + "'\t" + str(kapp_minimal_assumptions[rxn]))
        elif kapp_minimal_assumptions[rxn] <= kapp_ma_cutoff:
            kapp_ma_text.append("'" + rxn + "'\t" + str(kapp_ma_default))
        kapp_txt.append("'" + rxn + "'\t" + str(kapps_rba[rxn]))
with open('./kapp_default_values.txt', 'w') as f:
    f.write('\n'.join(output_info))
with open('./kapps_RBA.txt', 'w') as output:
    output.write('/\n' + '\n'.join(kapp_txt) + '\n' + '/')
with open('./kapps_minimal_assumptions.txt', 'w') as output:
    output.write('/' + '\n' + '\n'.join(kapp_ma_text) + '\n' + '/')
# if any rxns are inactive, we've predicted their kapps and know their fluxes
# use these to infer ENZSYN fluxes (ENZSYN=mu*flux/kapp)
kapp_test_text = []
if len(rxns_essential_NP) > 0:
    # # load fluxes
    # with open('min_flux_violation/min_flux_violation.flux.txt') as f:
    #     fluxes = f.read().split('\n')
    # fluxes = [i.split('\t') for i in fluxes]
    # fluxes = {i[0]:float(i[1]) for i in fluxes}
    # # calculate ENZSYN fluxes
    # enzsyn_fluxes = []
    c = 0
    for rxn in rxns_essential_NP:
        if rxn in kapps_rba.keys():
            # assume kapp is at its max value, to prevent excessive constraints
            kapps_rba[rxn] = kapp_max 
            kapp_test_text.append('Equation EnzCap'+str(c)+'; EnzCap'+str(c)+".. v('"+rxn.replace('RXN-','ENZLOAD-')+"')*"+str(kapps_rba[rxn])+" =e= %mu% * v('"+rxn+"');")
            c += 1
        # add its kapp as a constraint on ENZLOAD and flux
    #     if rxn in kapps_rba.keys() and rxn in fluxes.keys():
    #         enzsyn_fluxes.append('ENZSYN-' + rxn + '\t' + str(mu*fluxes[rxn]/kapps_rba[rxn]))
    # with open('./enzsyn_fluxes.txt', 'w') as f:
    #     f.write('\n'.join(enzsyn_fluxes))
    with open('./kapp_test.txt', 'w') as f:
        f.write('\n'.join(kapp_test_text))

# test if kapps work
#run_setting_file_from = './GAMS_setting_files/test_kapp_GAMS_settings.txt'
#run_setting_file_to = './min_flux_violation/test_kapp_GAMS_settings.txt'
#shutil.copy(run_setting_file_from, run_setting_file_to);
# shutil.copyfile('../../../../GAMS/runRBA.gms', './runRBA.gms')
# run RBA
os.system('cd min_flux_violation; module load gams; gams test_kapp.gms --nscale=' + str(nscale2) + output_redirect_str)
# check if RBA ran successfully
stop_if_run_failed('./min_flux_violation/test_kapp.modelStat.txt')
