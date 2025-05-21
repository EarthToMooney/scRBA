# parsimonious Protein FBA (pPFBA) method: minimize protein synthesis flux,
## recalculating kapps until they change negigibly
# update model-specific settings in kapp_options.py
from kapp_options import *
min_flux_cutoff = 0
# tolerance value added to assist with rounding kapps; set to 0 if not needed
tol = 0
# if kapp_diff_tol <= abs(1 - (new_kapp/old_kapp)), difference is considered significant
kapp_diff_tol = 0.9
# remove pPFBA_kapps_RBA.txt if it exists, to reset the files 
if True: 
    for file in ['./pPFBA_kapps_per_hr_without_unused_rxns.txt','./pPFBA_kapps_RBA.txt']:
        if os.path.exists(file):
            os.remove(file)
        # copy the pre-pPFBA counterparts of the pPFBA files
        shutil.copy(file.replace('pPFBA_',''),file)

attempt=0
significant_kapp_change = True
# Load enzyme info
df_enz = pd.read_excel('../../../../build_model/model/ENZYME_stoich_curation.tsv')

# Load path of set4 (enzyme-reaction many-to-many mapping)
set4_path = '../kapp_ambiguousLoad_case_resolve_common.txt'
# set4 handles cases of multiple enzymes corresponding to multiple rxns (many-to-many)
# columns: 
#   enzyme (for all enz that correspond to a rxn)
#   rxn: 
#       disableManually for any I don't want to count
#       noRxnFlux_resolveByManualCheck: behaves identically to disableManually, but perhaps serves as a warning
#   rule: 
#       divmaxenz for cases where I want to treat the [enzyme] as equaling its maximum value
gamsscaled_flux_path = './min_flux_violation/pPFBA.flux_gamsscaled.txt'
flux_path = './min_flux_violation/pPFBA.flux.txt'

while significant_kapp_change:
    significant_kapp_change = False # must find such a change to keep the loop going
    attempt+=1
    print('pPFBA round',attempt)
    os.system('cd min_flux_violation; module load gams; gams pPFBA.gms --nscale=' + str(nscale) + output_redirect_str)
    stop_if_run_failed('./min_flux_violation/pPFBA.modelStat.txt')

    res_metab = RBA_result(biom_id=biom_id)
    # res_metab.load_raw_flux('./min_flux_violation/pPFBA.flux_unscaled.txt')
    # extract lines starting with 'ENZSYN from fluxes
    with open(gamsscaled_flux_path) as f:
        rxns_used = f.read().split('\n')
    # find flux of total protein going into biomass: BIOSYN-PROTTOBIO
    # make dict of all PROWASTE rxns (except PROWASTE-TOTALPROTEIN) with name as key and flux as value
    prowaste_used = dict()
    prosyn_used = dict()
    enzsyn_used = []

    with open(flux_path,'w') as f:
        for i in rxns_used:
            i_twocol = i.replace('	v	 ','\t').replace('\tv\t','\t')
            if i_twocol == '':
                continue
            # print(i_twocol)
            rxn,v = i_twocol.split('\t')
            v = float(v)/nscale
            f.write('\t'.join([rxn,'v',str(v)])+'\n')
            if i_twocol[:8] == 'PROWASTE' and i_twocol.split('\t')[0] != 'PROWASTE-TOTALPROTEIN':
                prowaste_used[rxn] = v
            elif i_twocol[:7] == 'PROSYN-':
                prosyn_used[rxn] = v
            elif i_twocol[:7] == 'ENZSYN-':
                enzsyn_used.append(rxn+'\t'+str(v))
            elif i_twocol.split('\t')[0] == 'BIOSYN-PROTTOBIO':
                v_prot_to_bio = v

    # Add new set of constraints to reflect unused but observed protein production, similar to dummy protein (for every prot, v('PROSYN-prot') >= prowaste_used('PROWASTE-prot') / v_prot_to_bio * %nscale%)

    ## UNUSED: make separate file for defaults, so not every options file needs to be updated when the model changes
    ## pro_measured_unused_default_path = '../../../../prosyn_constraints_measured_but_unused_default.txt'
    ## if os.path.exists(pro_measured_unused_default_path):

    # pro_measured_unused_path = '../../../../prosyn_constraints_measured_but_unused_'+dir_name+'.txt'
    # c = 0
    # with open(pro_measured_unused_path,'w') as f:
    #     # for each protein in PROWASTE-used, allow flux through PROWASTE-PROT to be as high as infinity (to avoid scaling issues)
    #     for prot in prowaste_used.keys():
    #         prot_id = prot.replace('PROWASTE-','')
    #         prosyn_id = 'PROSYN-'+prot_id
    #         if prosyn_id in prosyn_used.keys():
    #             f.write("v.up('"+prosyn_id+"') = inf; v.up('"+prot+"') = inf;\nEquation minProt"+str(c)+"; minProt"+str(c)+".. v('"+prosyn_id+"') =g= v('BIOSYN-PROTTOBIO') * "+str(prowaste_used[prot])+" / "+str(v_prot_to_bio)+";\n")
    #             c += 1
    # write to new file
    kapp_input_enzsyn_path = './pPFBA_kapp_input_enzsyn.txt'
    with open(kapp_input_enzsyn_path, 'w') as f:
        f.write('\n'.join(enzsyn_used))
    res_metab.load_raw_flux(flux_path)
    res_metab.calculate_metabolic_flux()

    res_esyn = RBA_result(biom_id=biom_id, twocol_format=True)
    # load_raw_flux using a version of min_flux_violation.flux.txt that only contains lines starting with 'ENZSYN
    #res_esyn.load_raw_flux('./min_flux_violation/enz_alloc.flux.txt')
    # res_esyn.load_raw_flux('./min_flux_violation/enz_flux_calculation.txt')
    res_esyn.load_raw_flux(kapp_input_enzsyn_path)
    with open("../../../../GAMS/model/RBA_rxns_rxnmetabolicnetwork.txt","r") as f: 
        all_rxns = [line.strip().replace("'","") for line in f.readlines() if line.strip() != "/"]
    with open("../../../../GAMS/model/GSM_rxn_ids.txt","r") as f: 
        all_GSM_rxns = [line.strip().replace("'","") for line in f.readlines() if line.strip() != "/"]

    mu = res_metab.growth_rate
    print('Growth rate:', mu)
    print('EX_glc__D_e', res_metab.metabolic_flux['EX_glc__D_e'])
    print('EX_o2_e', res_metab.metabolic_flux['EX_o2_e'])
    for k,v in res_metab.metabolic_flux.items():
        if k[:3] == 'EX_' and v > 1e-6:
            print(k,v)
    #### Map rxn to enz
    rxndict = {k:[] for k,(v) in res_metab.metabolic_flux.items() if abs(v) > 0}

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
    rxndict_zeroCost = {k:v for k,v in rxndict.items() if v == {'zeroCost'}}
    rxndict = {k:v for k,v in rxndict.items() if v != {'zeroCost'} and 'zeroCost' not in v}
    enzdict = {k:set(v) for k,v in enzdict.items() if 'zeroCost' not in [k,v] and 'zeroCost' not in v}
    # Find enzymes that together carry a total load of reactions
    x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
    enz_share_rxn_load = set().union(*[v for v in x.values() if v != 'zeroCost'])

    # Find enzymes that individually carry loads of multiple reactions
    x = {k:v for k,v in enzdict.items() if len(v) > 1.5}
    enz_multiload = set(x.keys())

    # Set 1: Enzyme-reaction one-to-one load mapping
    set1 = set(enzdict.keys()) - enz_share_rxn_load - enz_multiload
    set1 = set([i for i in set1 if len(enzdict[i]) > 0.5])
    # write set1 to file
    with open('./pPFBA_kapp_set1_one_to_one.txt', 'w') as f:
        f.write('\n'.join(set1))

    # Set 2: Enzyme-reaction one-to-many
    set2 = enz_multiload - enz_share_rxn_load
    # write set2 to file
    with open('./pPFBA_kapp_set2_one_to_many.txt', 'w') as f:
        f.write('\n'.join(set2))

    # Set 3: Enzyme-reaction many-to-one
    set3 = enz_share_rxn_load - enz_multiload
    # write set3 to file
    with open('./pPFBA_kapp_set3_many_to_one.txt', 'w') as f:
        f.write('\n'.join(set3))

    # Set 4: Enzyme-reaction many-to-many
    set4 = enz_share_rxn_load & enz_multiload
    # write set4 to file
    with open('./pPFBA_kapp_set4_many_to_many.txt', 'w') as f:
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
        if i not in rxns_to_ignore_for_kapps:
            print(i)
    # Find enzymes whose together carry a total load of reactions
    x = {k:v for k,v in rxndict.items() if len(v) > 1.5}
    enz_share_rxn_load = set().union(*[v for v in x.values() if v != 'zeroCost'])

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
        
    # From manual checking of the sets, some enzymes in set2 and set3 require special treatment in calculation
    # which are recorded in set4. Thus excluded those in set2 and set3
    set2 = set2 - set4 - {'zeroCost'}
    set3 = set3 - set4 - {'zeroCost'}

    kapp = dict()
    kapps_rba = dict()
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
    with open('../exclude_parameterization_list.txt') as f:
        excl = f.read().split('\n')
    # # TEST: add all rxns with fluxes below min_flux_cutoff to excl
    # for k,v in res_metab.metabolic_flux.items():
    #    if abs(v) < min_flux_cutoff:
    #        excl.append(k)
    excl = [r for r in excl if r != '']
    excl += []

    texts = ['rxnid\tkapp (1/s)']
    perhr_texts_used_only = ['/']
    for k,v in kapp.items():
        if k not in excl:
            texts.append(k + '\t' + str(v))
        
    with open('./pPFBA_kapps_in_vivo.txt', 'w') as f:
        f.write('\n'.join(texts))
    kapp_med = np.median(list(kapp.values())) * 3600
    kapp_max = max(list(kapps_rba.values())) 
    kapp_ma_med = np.median(list(kapp_minimal_assumptions.values()))
    try:
        kapp_ma_max = np.max(list(kapp_minimal_assumptions.values()))
        kapp_ma_default = kapp_ma_med
    except:
        kapp_ma_max = 0
        kapp_ma_default = kapp_med
    default_kapp = kapp_med + tol
    # return kapps in a format suitable for use in RBA
    kapp_txt = ['/']

    # check if any inactive rxns were used in B2
    with open('min_flux_violation/enz_alloc.flux_essential_inactive_rxns.txt') as f:
        rxns_essential_inactive = f.read().split('\n')
    rxns_essential_inactive = [i for i in rxns_essential_inactive if i != '']
    rxns_essential_inactive = [i.split('\t')[0] for i in rxns_essential_inactive]
    for r in rxns_essential_inactive:
        if r not in kapps_rba.keys():
            kapps_rba[r] = default_kapp + tol
    for k,v in kapps_rba.items():
        perhr_texts_used_only.append("'" + k + "'" + '\t' + str(v))
    path='./pPFBA_kapps_per_hr_without_unused_rxns.txt'
    with open(path, 'r') as f:
        old_kapps = dict()
        new_kapps = dict()
        # read old version of file (to filter out non-default and irrelevant values) and compare to new version
        # if file isn't empty
        if os.path.getsize(path) > 0:
            #for dic, txt in {old_kapps: f.read().split('\n'), new_kapps: perhr_texts_used_only}.items():
            for i,txt in enumerate([f.read().split('\n'), perhr_texts_used_only]):
                #print('dic:',str(dic),txt)
                for line in txt:
                    if '/' not in line.split('\t'):
                        k, v = line.split('\t')
                        if i == 0:
                            old_kapps[k.replace("'","")] = v
                        else:
                            new_kapps[k.replace("'","")] = v
            #print('old kapps:',str(old_kapps))
            # print any differences
            for k,v in old_kapps.items():
                if k in new_kapps.keys():
                    if kapp_diff_tol <= abs(1-(float(new_kapps[k])/float(v))):
                        print('kapp changed by at least',kapp_diff_tol*100,"% :",k,'old:',v,'new:',new_kapps[k])
                        significant_kapp_change = True
                else:
                    print('kapp removed:',k,v)
                    significant_kapp_change = True
        # old_kapps = {k:v for k,v in f.read().split('\n').split('\t') if '/' not in [k,v]}
        # new_kapps = {k:v for k,v in perhr_texts_used_only.split('\t') if '/' not in [k,v]}
    with open('./pPFBA_kapps_per_hr_without_unused_rxns_round'+str(attempt)+'.txt', 'w') as f:
        f.write('\n'.join(sorted(perhr_texts_used_only) + ['/']))
    with open(path, 'w') as f:
        f.write('\n'.join(sorted(perhr_texts_used_only) + ['/']))

    output_info = []
    # add all other rxns with default kapp
    for rxn in all_rxns: 
        enz = rxn.split('-')[-1]
        rxn_name_without_enz = rxn.split('-')[-2]
        if enz not in spont_rxn_suffixes:
            if rxn in kapps_rba.keys():
                if kapps_rba[rxn] <= min_flux_cutoff:
                    output_info.append('kapp <= cutoff for ' + rxn)
                    kapps_rba[rxn] = default_kapp
            elif rxn in rxns_essential_inactive:
                kapps_rba[rxn] = default_kapp # to minimize the risk of kapps making growth infeasible 
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
    with open('./pPFBA_kapp_default_values.txt', 'w') as f:
        f.write('\n'.join(output_info))
    with open('./pPFBA_kapps_RBA.txt', 'w') as output:
        output.write('\n'.join(kapp_txt) + '\n' + '/')
    with open('./pPFBA_kapps_minimal_assumptions.txt', 'w') as output:
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
        with open('./pPFBA_kapp_test.txt', 'w') as f:
            f.write('\n'.join(kapp_test_text))

    # test if kapps work
    from simulate import get_GAMS_modelStat
    #run_setting_file_from = './GAMS_setting_files/test_kapp_GAMS_settings.txt'
    #run_setting_file_to = './min_flux_violation/test_kapp_GAMS_settings.txt'
    #shutil.copy(run_setting_file_from, run_setting_file_to);
    # shutil.copyfile('../../../../GAMS/runRBA.gms', './runRBA.gms')
    # run RBA
    os.system("cd min_flux_violation; module load gams; gams test_kapp.gms --kapp_path='../pPFBA_kapps_RBA.txt' --nscale=" + str(nscale2) + output_redirect_str)
    # check if RBA ran successfully
    stop_if_run_failed('./min_flux_violation/test_kapp.modelStat.txt')

