import cobra
# based on initial file from Hoang Viet Dinh
# updated to Python 3

def test_import():
    print('Import ok, common custom functions')

def compile_dictionary_from_text(fpath, sep='\t', keypos=0, valuepos=1, skiprows=0):
    """Each dictionary key:value pair is assumed to be on a line. Pairs was separated by newline (backslash n).
    fpath: directory to text file
    sep: delimiter used in text file
    keypos: position (start from 0) for dictionary key
    valuepos: position (start from 0) for dictionary value"""
    from collections import OrderedDict

    with open(fpath) as f:
        x = f.read().split('\n')[skiprows:]
    outdict = OrderedDict()
    for i in x:
        if i != '':
            i2 = i.split(sep)
            outdict[i2[keypos]] = i2[valuepos]
    return outdict

def make_escher_csv(mflux, path):
    import csv
    import pandas
    if isinstance(mflux, pandas.core.series.Series) != True:
        mflux = mflux.fluxes
    with open(path, 'w') as f:
        fcsv = csv.writer(f, delimiter=',')
        fcsv.writerow(['Rxn', 'Flux'])
        for rxn in mflux.index:
            fcsv.writerow([rxn, mflux[rxn]])

def duplicate_metabolite(metabolite, modelfrom, modelto):
    """Duplicate a 'metabolite' from 'modelfrom' to 'modelto'"""
    if isinstance(metabolite, str):
        metfrom = modelfrom.metabolites.get_by_id(metabolite)
    else:
        metfrom = metabolite
        
    import cobra
    metto = cobra.Metabolite(metfrom.id)
    modelto.add_metabolites([metto])
    
    metto = modelto.metabolites.get_by_id(metfrom.id)
    metto.name = metfrom.name
    metto.formula = metfrom.formula
    metto.charge = metfrom.charge
    metto.compartment = metfrom.compartment
    metto.annotation = metfrom.annotation
    metto.notes = metfrom.notes
    
    return None

def duplicate_reaction(reaction, modelfrom, modelto):
    """Duplicate a 'reaction' from 'modelfrom' to 'modelto'"""
    if isinstance(reaction, str):
        rxnfrom = modelfrom.reactions.get_by_id(reaction)
    else:
        rxnfrom = reaction
        
    mets_in_modelto = [met.id for met in modelto.metabolites]
    for met in rxnfrom.metabolites:
        if met.id not in mets_in_modelto:
            duplicate_metabolite(met, modelfrom, modelto)
    
    import cobra
    rxnto = cobra.Reaction(rxnfrom.id)
    modelto.add_reactions([rxnto])
    
    rxnto = modelto.reactions.get_by_id(rxnfrom.id)
    rxnto.name = rxnfrom.name
    rxnto.reaction = rxnfrom.reaction
    rxnto.bounds = rxnfrom.bounds
    rxnto.gene_reaction_rule = rxnfrom.gene_reaction_rule
    rxnto.annotation = rxnfrom.annotation
    rxnto.notes = rxnfrom.notes
    
    return None
    
def transfer_gene_info(gene, modelfrom, modelto):
    """Duplicate a 'gene' from 'modelfrom' to 'modelto'"""
    if isinstance(gene, str):
        genefrom = modelfrom.genes.get_by_id(gene)
    else:
        genefrom = gene
        
    import cobra
    genes_in_modelto = [g.id for g in modelto.genes]
    
    if gene not in genes_in_modelto:
        print(gene + ' is not in the modelto')
    else:
        geneto = modelto.genes.get_by_id(gene)
        geneto.name = genefrom.name
        geneto.annotation = genefrom.annotation
        
    return None

def query_metabolites(model, keyword, category):
    """Search for 'keyword' in specific 'category' in 'model'"""
    output = []
    if category in ['KEGG', 'kegg', 'kegg.compound']:
        for met in model.metabolites:
            if 'kegg.compound' in met.annotation.keys() and met.annotation['kegg.compound'] == keyword:
                output.append(met.id)
                print(met.id, met.name)
                
    elif category in ['name']:
        for met in model.metabolites:
            if keyword.lower() in met.name.lower():
                output.append(met.id)
                print(met.id, met.name)
                
    elif category in ['id', 'ID']:
        for met in model.metabolites:
            if keyword.lower() in met.id.lower():
                output.append(met.id)
                print(met.id, met.name)
                
    return output

def query_reactions(model, keyword, category):
    """Search for 'keyword' in specific 'category' in 'model'"""
    output = []
    if category in ['KEGG', 'kegg', 'kegg.reaction']:
        for rxn in model.reactions:
            if 'kegg.reaction' in rxn.annotation.keys() and rxn.annotation['kegg.reaction'] == keyword:
                output.append(rxn.id)
                print(rxn.id, rxn.name)
                
    elif category in ['EC', 'ec' , 'EC number', 'ec number', 'ecnum', 'ec-code']:
        for rxn in model.reactions:
            if 'ec-code' in rxn.annotation.keys() and rxn.annotation['ec-code'] == keyword:
                output.append(rxn.id)
                print(rxn.id, rxn.name)
                
    elif category in ['name']:
        for rxn in model.reactions:
            if keyword.lower() in rxn.name.lower():
                output.append(rxn.id)
                print(rxn.id, rxn.name)
                
    elif category in ['id', 'ID']:
        for rxn in model.metabolites:
            if keyword.lower() in rxn.id.lower():
                output.append(rxn.id)
                print(rxn.id, rxn.name)
    return output

def generate_metabolite_in_diff_compartment(modelto, metid, compto, modelfrom, met_copy_source, verbose=False):
    """Generate new metabolite ('metid') in 'modelto' with the same activity with 'rxn_copy_source'
    from 'modelfrom' but in different compartment, 'compto'."""
    
    mets_model = [r.id for r in modelto.metabolites]
    if metid in mets_model:
        print(metid + ' (a metabolite) already in the model')
        
    else:
        metfrom = modelfrom.metabolites.get_by_id(met_copy_source)

        met = cobra.Metabolite(metid)
        modelto.add_metabolites([met])
        met = modelto.metabolites.get_by_id(metid)
        met.name = str(metfrom.name)
        met.formula = str(metfrom.formula)
        met.charge = float(metfrom.charge)
        met.compartment = compto
        
        notrans = ['yeast_id', 'yeastkegg_id'] # Unique for my case (July 07 2018, should be removed for code clarity, but otherwise safe to just let it be)
        anns = {k:v for k,v in metfrom.annotation.items() if k not in notrans}
        met.annotation = anns
        
        if verbose:
            print(metid + ' (a metabolite) is added to the model')
    
    return modelto

def generate_reaction_in_diff_compartment(modelto, rxnid, compto, modelfrom, rxn_copy_source, verbose=False):
    """Generate new reaction ('rxnid') in 'modelto' with the same activity with 'rxn_copy_source'
    from 'modelfrom' but in different compartment, 'compto'. WARNING: this function will convert
    all original metabolites' compartments to a single 'compto' compartment"""
    
    rxns_model = [r.id for r in modelto.reactions]
    if rxnid in rxns_model:
        print(rxnid + ' (a reaction) already in the model')
    else:
        rxnfrom = modelfrom.reactions.get_by_id(rxn_copy_source)

        rxn = cobra.Reaction(rxnid)
        modelto.add_reactions([rxn])
        rxn = modelto.reactions.get_by_id(rxnid)
        rxn.name = str(rxnfrom.name)
        rxn.gene_reaction_rule = rxnfrom.gene_reaction_rule

        stoich = dict()
        for metfrom, coeff in rxnfrom.metabolites.items():
            compmet = metfrom.id.split('_')[-1]
            idcore = metfrom.id[:-len(compmet) - 1]
            metid = idcore + '_' + compto
            stoich[metid] = coeff

            mets_model = [m.id for m in modelto.metabolites]
            if metid not in mets_model:
                modelto = generate_metabolite_in_diff_compartment(modelto, metid, compto, modelfrom, metfrom.id)

        stoich_new = {modelto.metabolites.get_by_id(k):v for k,v in stoich.items()}
        rxn.add_metabolites(stoich_new)
        
        notrans = ['yeast_id', 'yeastkegg_id'] # Unique for my case (July 07 2018, should be removed for code clarity, but otherwise safe to just let it be)
        anns = {k:v for k,v in rxnfrom.annotation.items() if k not in notrans}
        rxn.annotation = anns

        if verbose:
            print('\t'.join([rxn.id, rxn.reaction]))
    
    return modelto

def execute_command(model, model_donor, df_cmds, verbose=False):
    attr_dict = {i:i for i in ['id', 'reaction', 'name', 'subsystem', 'lower_bound',
                               'upper_bound', 'compartment', 'notes', 'formula', 'charge']}
    attr_dict['gpr'] = 'gene_reaction_rule'
    
    comps_model = list(model.compartments.keys())
    rxns_model_donor = [rxn.id for rxn in model_donor.reactions]
    mets_model_donor = [met.id for met in model_donor.metabolites]
    
    for i in df_cmds.index:
        rxns = model.reactions
        mets = model.metabolites

        cmd = df_cmds.command[i]
        objid = df_cmds.id[i]
        obj = df_cmds.object_type[i]

        rxns_model = [rxn.id for rxn in model.reactions]
        mets_model = [met.id for met in model.metabolites]
        genes_model = [g.id for g in model.genes]
        
        if obj.lower() in ['r', 'rxn', 'reaction', 'reactions']:
            if cmd in ['create', 'make', 'generate']:
                if objid in rxns_model:
                    print('\t'.join([objid, obj, 'already exists in the model']))
                else:
                    model.add_reactions([cobra.Reaction(objid)])
            elif cmd == 'retrieve':
                duplicate_reaction(objid, model_donor, model)
            elif cmd == 'remove':
                model.remove_reactions([objid])
            else:
                cmd1 = cmd.split(':')[0]
                cmd2 = cmd[len(cmd1)+1:]
                if cmd1 == 'copy':
                    compto = objid.split('_')[-1]
                    if compto not in comps_model:
                        compto = 'c'
                    if cmd2 in rxns_model:
                        model = generate_reaction_in_diff_compartment(model, objid, compto, model, cmd2)
                    elif cmd2 in rxns_model_donor:
                        model = generate_reaction_in_diff_compartment(model, objid, compto, model_donor, cmd2)
                    else:
                        print('When try to copy ' + objid + ' in different compartment, ' + cmd2 + ' is not in either model or model_donor')

                elif cmd1 in ['id', 'reaction', 'name', 'subsystem', 'gpr', 'lower_bound', 'upper_bound']:
                    if objid not in rxns_model:
                        print(objid + ' is not in the model as a reaction, cannot change ' + cmd1)
                    else:
                        if cmd1 in ['id', 'reaction', 'name', 'subsystem', 'gpr']:
                            setattr(model.reactions.get_by_id(objid), attr_dict[cmd1], cmd2)
                        elif cmd1 in ['lower_bound', 'upper_bound']:
                            setattr(model.reactions.get_by_id(objid), attr_dict[cmd1], float(cmd2))
                        
                elif cmd1 == 'annotation':
                    if objid not in rxns_model:
                        print(objid + ' is not in the model as a reaction, cannot change ' + cmd1)
                    else:
                        if cmd2 == 'remove':
                            model.reactions.get_by_id(objid).annotation = dict()
                        else:
                            if cmd2.split('-')[0] == 'remove_part':
                                if cmd2.split('-')[1] in model.reactions.get_by_id(objid).annotation.keys():
                                    del model.reactions.get_by_id(objid).annotation[cmd2.split('-')[1]]
                            elif cmd2.split('-')[0] in ['ec', 'EC']:
                                model.reactions.get_by_id(objid).annotation['ec-code'] = cmd2.split('-')[1]
                            else:
                                cmd2_part1 = cmd2.split('-')[0]
                                model.reactions.get_by_id(objid).annotation[cmd2_part1] = cmd2[len(cmd2_part1)+1:]
                
                elif cmd1 == 'notes':
                    if objid not in rxns_model:
                        print(objid + ' is not in the model as a reaction, cannot change ' + cmd1)
                    else:
                        if cmd2 == 'remove':
                            model.reactions.get_by_id(objid).notes = dict()
                        else:
                            if cmd2.split('-')[0] == 'remove_part':
                                if cmd2.split('-')[1] in model.reactions.get_by_id(objid).notes.keys():
                                    del model.reactions.get_by_id(objid).notes[cmd2.split('-')[1]]
                            elif cmd2.split('-')[0] == 'notes':
                                if 'notes' in model.reactions.get_by_id(objid).notes.keys():
                                    note = model.reactions.get_by_id(objid).notes['notes']
                                    model.reactions.get_by_id(objid).notes['notes'] = '. '.join([note, str(cmd2.split('-')[1])])
                                else:
                                    model.reactions.get_by_id(objid).notes['notes'] = str(cmd2.split('-')[1])
                            else:
                                model.reactions.get_by_id(objid).notes[cmd2.split('-')[0]] = cmd2.split('-')[1]
        
        if obj.lower() in ['m', 'met', 'metabolite', 'metabolites', 'compound']:
            if cmd in ['create', 'make', 'generate']:
                if objid in mets_model:
                    print('\t'.join([objid, obj, 'already exists in the model']))
                else:
                    model.add_metabolites([cobra.Metabolite(objid)])
            if cmd == 'retrieve':
                duplicate_metabolite(objid, model_donor, model)
            elif cmd == 'remove':
                model.remove_metabolites([model.metabolites.get_by_id(objid)])
            else:
                cmd1 = cmd.split(':')[0]
                cmd2 = cmd[len(cmd1)+1:]
                if cmd1 == 'copy':
                    compto = objid.split('_')[-1]
                    if cmd2 in mets_model:
                        model = generate_metabolite_in_diff_compartment(model, objid, compto, model, cmd2)
                    elif cmd2 in mets_model_donor:
                        model = generate_metabolite_in_diff_compartment(model, objid, compto, model_donor, cmd2)
                    else:
                        print('When try to copy ' + objid + ' in different compartment, ' + cmd2 + ' is not in either model or model_donor')

                elif cmd1 in ['id', 'name', 'compartment', 'formula', 'charge']:
                    if objid not in mets_model:
                        print(objid + ' is not in the model as a metabolite, cannot change ' + cmd1)
                    else:
                        if cmd1 in ['id', 'name', 'compartment', 'formula']:
                            setattr(model.metabolites.get_by_id(objid), attr_dict[cmd1], cmd2)
                        elif cmd1 in ['charge']:
                            setattr(model.metabolites.get_by_id(objid), attr_dict[cmd1], float(cmd2))
                        
                elif cmd1 == 'annotation':
                    if objid not in mets_model:
                        print(objid + ' is not in the model as a metabolite, cannot change ' + cmd1)
                    else:
                        if cmd2 == 'remove':
                            model.metabolites.get_by_id(objid).annotation = dict()
                        else:
                            if cmd2.split('-')[0] == 'remove_part':
                                if cmd2.split('-')[1] in model.metabolites.get_by_id(objid).annotation.keys():
                                    del model.metabolites.get_by_id(objid).annotation[cmd2.split('-')[1]]
                            else:
                                cmd2_part1 = cmd2.split('-')[0]
                                model.metabolites.get_by_id(objid).annotation[cmd2_part1] = cmd2[len(cmd2_part1)+1:]
                                
                elif cmd1 == 'notes':
                    if objid not in mets_model:
                        print(objid + ' is not in the model as a reaction, cannot change ' + cmd1)
                    else:
                        if cmd2 == 'remove':
                            model.metabolites.get_by_id(objid).notes = dict()
                        else:
                            if cmd2.split('-')[0] == 'remove_part':
                                if cmd2.split('-')[1] in model.metabolites.get_by_id(objid).notes.keys():
                                    del model.metabolites.get_by_id(objid).notes[cmd2.split('-')[1]]
                            elif cmd2.split('-')[0] == 'notes':
                                if 'notes' in model.metabolites.get_by_id(objid).notes.keys():
                                    note = model.metabolites.get_by_id(objid).notes['notes']
                                    model.metabolites.get_by_id(objid).notes['notes'] = '. '.join([note, str(cmd2.split('-')[1])])
                                else:
                                    model.metabolites.get_by_id(objid).notes['notes'] = str(cmd2.split('-')[1])
                            else:
                                model.metabolites.get_by_id(objid).notes[cmd2.split('-')[0]] = cmd2.split('-')[1]
        
        if obj.lower() in ['g', 'gene', 'genes']:
            if cmd == 'retrieve':
                transfer_gene_info(objid, model_donor, model)
            else:
                cmd1 = cmd.split(':')[0]
                cmd2 = cmd[len(cmd1)+1:]
                if cmd1 in ['id', 'name']:
                    if objid not in genes_model:
                        print(objid + ' is not in the model as a gene, cannot change ' + cmd1)
                    else:
                        if cmd1 == 'id':
                            cobra.manipulation.modify.rename_genes(model, {objid:cmd2})
                        else:
                            setattr(model.genes.get_by_id(objid), attr_dict[cmd1], cmd2)

                elif cmd1 == 'annotation':
                    if objid not in genes_model:
                        print(objid + ' is not in the model as a gene, cannot change ' + cmd1)
                    else:
                        if cmd2 == 'remove':
                            model.genes.get_by_id(objid).annotation = dict()
                        else:
                            if cmd2.split('-')[0] == 'remove_part':
                                if cmd2.split('-')[1] in model.genes.get_by_id(objid).annotation.keys():
                                    del model.genes.get_by_id(objid).annotation[cmd2.split('-')[1]]
                            else:
                                cmd2_part1 = cmd2.split('-')[0]
                                model.genes.get_by_id(objid).annotation[cmd2_part1] = cmd2[len(cmd2_part1)+1:]
                                
                elif cmd1 == 'notes':
                    if objid not in genes_model:
                        print(objid + ' is not in the model as a reaction, cannot change ' + cmd1)
                    else:
                        if cmd2 == 'remove':
                            model.genes.get_by_id(objid).notes = dict()
                        else:
                            if cmd2.split('-')[0] == 'remove_part':
                                if cmd2.split('-')[1] in model.genes.get_by_id(objid).notes.keys():
                                    del model.genes.get_by_id(objid).notes[cmd2.split('-')[1]]
                                elif cmd2.split('-')[0] == 'notes':
                                    if 'notes' in model.genes.get_by_id(objid).notes.keys():
                                        note = model.genes.get_by_id(objid).notes['notes']
                                        model.genes.get_by_id(objid).notes['notes'] = '. '.join([note, str(cmd2.split('-')[1])])
                                    else:
                                        model.genes.get_by_id(objid).notes['notes'] = str(cmd2.split('-')[1])
                            else:
                                model.genes.get_by_id(objid).notes[cmd2.split('-')[0]] = cmd2.split('-')[1]
        
        if verbose:
            print('\t'.join(df_cmds.loc[i, :]))
            
    return model

def compile_elements_from_formula(formula):
    import re
    from collections import OrderedDict
    elems = re.findall('[A-Z][a-z]*', formula)
    elem_dict = OrderedDict()
    for elem in elems:
        # handles positive and negative coefficients (negative for adjust_formula function)
        elemval = re.search('(?<=' + elem + ')[-\d.]*', formula).group()
        if elemval == '':
            elem_dict[elem] = 1.0
        else:
            elem_dict[elem] = float(elemval)
    return elem_dict

def compile_formula_from_elements(elem_dict):
    from collections import OrderedDict
    formula = ''
    for k,v in elem_dict.items():
        if v == 0:
            continue
        elif v == 1:
            formula += k
        elif type(v) == int or v.is_integer():
            formula += k + str(int(v))
        else:
            formula += k + str(v)
    return formula

def adjust_formula(formula_to_change, formula_to_change_by):
    '''Add (or subtract, by using negative coefficients) atoms from formula_to_change_by to formula_to_change.
    \nIt helps to put longer elements later (e.g., Co after C) since the results may be incorrect otherwise.
    '''
    import cobra
    elem_dict1 = compile_elements_from_formula(formula_to_change)
    elem_dict2 = compile_elements_from_formula(formula_to_change_by)
    elem_dict3 = dict()
    for k,v in elem_dict1.items():
        if k in elem_dict2.keys():
            newcoeff = v + elem_dict2[k]
            if newcoeff > 0:
                elem_dict3[k] = v + elem_dict2[k]
            elif newcoeff == 0:
                continue
            else:
                print(formula_to_change + ' + ' + formula_to_change_by + ' yields negative coefficient for element ' + k + '; removing it instead')
        else:
            elem_dict3[k] = v
    for k,v in elem_dict2.items():
        if k not in elem_dict1.keys():
            elem_dict3[k] = v
    return compile_formula_from_elements(elem_dict3)

def check_mass_balance_cobra(reaction, model):
    import cobra, re
    from collections import OrderedDict
    
    if type(reaction) in [str, str]:
        rxn = model.reactions.get_by_id(reaction)
    else:
        rxn = reaction
    
    elem_mets = {met.id:compile_elements_from_formula(met.formula) for met in rxn.metabolites.keys()}
    coeffs = {met.id:coeff for met,coeff in rxn.metabolites.items()}
    imbal = OrderedDict()
    for v in elem_mets.values():
        for elem in v.keys():
            if elem not in elem_mets.keys():
                imbal[elem] = 0
    
    for metid, v in elem_mets.items():
        for elem, elemval in v.items():
            imbal[elem] += coeffs[metid] * elemval
            
    imbal['charge'] = 0
    for met in rxn.metabolites.keys():
        imbal['charge']  += coeffs[met.id] * met.charge
        
    return imbal

def tab_print_adjustment(list_in, char_length=0):
    clen = max([max([len(str(i)) for i in list_in]), char_length])
    list_out = [str(i) + ' '*(clen - len(str(i))) for i in list_in]
    return list_out

def make_cobra_model_from_excel(excelFile, sheetDict, propDict, modelName='model'):
    import pandas as pd
    import cobra
    
    dfMets = pd.read_excel(excelFile, sheet_name=sheetDict['metabolites'])
    dfMets = dfMets.rename({v:k for k,v in propDict['metabolites'].items()}, axis='columns')
    dfRxns = pd.read_excel(excelFile, sheet_name=sheetDict['reactions'])
    dfRxns = dfRxns.rename({v:k for k,v in propDict['reactions'].items()}, axis='columns')

    metPropsOptional = ['formula', 'charge']
    rxnPropsOptional = ['gene_reaction_rule', 'subsystem', 'reversibility', 'lower_bound',
                        'upper_bound', 'objective_coefficient']

    model = cobra.Model(modelName)
    for i in dfMets.index:
        model.add_metabolites([cobra.Metabolite(dfMets.id[i])])
        met = model.metabolites.get_by_id(dfMets.id[i])
        met.name = dfMets.name[i]
        for prop in metPropsOptional:
            if prop in dfMets.columns and pd.isnull(dfMets.loc[i, prop]) == False:
                setattr(met, prop, dfMets.loc[i, prop])

    for i in dfRxns.index:
        model.add_reactions([cobra.Reaction(dfRxns.id[i])])
        rxn = model.reactions.get_by_id(dfRxns.id[i])
        rxn.name = dfRxns.name[i]
        rxn.reaction = dfRxns.reaction[i]
        for prop in rxnPropsOptional:
            if prop in dfRxns.columns and pd.isnull(dfRxns.loc[i, prop]) == False:
                setattr(rxn, prop, dfRxns.loc[i, prop])
    
    return model

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def metabolites_dict_from_reaction_equation_RBA(eqn, split=False):
    import re
    rStr, pStr = re.split('-->|->|<--|<-|<=>|<->', eqn)
    rStr = rStr.strip(' ')
    pStr = pStr.strip(' ')
    rs = re.split(' \+ | \+|\+ |\+', rStr)
    ps = re.split(' \+ | \+|\+ |\+', pStr)

    r_dict = dict()
    for r in rs:
        if ' ' in r:
            val, met = re.split('\s+', r)
            r_dict[met] = -float(val)
        else:
            r_dict[r] = -1.0

    p_dict = dict()
    for p in ps:
        if ' ' in p:
            val, met = re.split('\s+', p)
            p_dict[met] = float(val)
        else:
            p_dict[p] = 1.0
            
    if split:
        return r_dict, p_dict
    else:
        return merge_two_dicts(r_dict, p_dict)

def noncomp_id(metid):
    comp = metid.split('_')[-1]
    return metid[:-len(comp)-1]

def build_reaction_equation_from_metabolites_dict_RBA(met_dict, arrow='<=>', floatdecimal=6):
    """Takes dict (keys=metabolite IDs, values=stoichiometric coefficients) and returns reaction equation string"""
    lhs = []; rhs = [];
    for k,v in met_dict.items():
        v = float(v)
        if v == -1:
            lhs.append(k)
        elif v == 1:
            rhs.append(k)
        elif v < 0 and v != -1 and v.is_integer():
            lhs.append(' '.join([str(-int(v)), k]))
        elif v > 0 and v != 1 and v.is_integer():
            rhs.append(' '.join([str(int(v)), k]))
        elif v < 0 and v != -1:
            lhs.append(' '.join([('{:.' + str(floatdecimal) + 'f}').format(-v), k]))
        elif v > 0 and v != 1:
            rhs.append(' '.join([('{:.' + str(floatdecimal) + 'f}').format(v), k]))
    return ' '.join([ ' + '.join(lhs), arrow, ' + '.join(rhs)])

def find_biomass_reactions(model, properties='id'):
    '''Finds biomass reactions in a model (SBO:0000629) and returns the properties you want (by default, ID).'''
    biomass_rxns = []
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower():
            biomass_rxns.append(rxn.properties)
            break
        elif 'sbo' in rxn.annotation:
            if 'SBO:0000629' in rxn.annotation['sbo']:
                biomass_rxns.append(rxn.properties)
                break
    return biomass_rxns

def find_pathways(model, reactant_id, product_id,exclude_common_mets=True,stop_at_first=True):
    """UNFINISHED: find pathways that allow a specified reactant to form a specified product in a model.
    \nBy default, to reduce time and improve accuracy, excludes water, lone protons, and metabolites commonly used as currency metabolites or cofactors (e.g. ATP, ADP, NADH, NADPH, etc.)."""
    met_from = model.metabolites.get_by_id(reactant_id)
    met_to = model.metabolites.get_by_id(product_id)
    if exclude_common_mets:
        excluded_met_names = ['H2O','H','H+','PROTON','ATP','ADP','AMP','NADH','NADPH','NAD+','NADP+','NAD','NADP','NADH2','NADPH2','NADP(+)','NADP(H)','NAD(H)','NAD(P)H','NAD(P)','FAD','FADH2','FMN','FMNH2']
        # find all mets with names matching those in the excluded list
        excluded_mets = [met for met in model.metabolites if met.name in excluded_met_names]
    pathways = []
    # find all products for met_to
    for rxn in met_to.reactions:
            # check if the rxn can produce met_to
            if rxn.upper_bound > 0 and met_to in rxn.products:
                if met_from in rxn.reactants:
                    pathways.append(rxn.id)
                    if stop_at_first:
                        break
                # else:
                    # add the reactants of the rxn to the list of mets to check
            elif rxn.lower_bound < 0 and met_to in rxn.reactants:
                if met_from in rxn.products:
                    pathways.append(rxn.id)
                    if stop_at_first:
                        break
    return pathways

def build_stoichiometry_string(obj, number_delim=':', metabolite_delim=','):
    if type(obj) == str:
        met_dict = metabolites_dict_from_reaction_equation_RBA(obj)
    if type(obj) == dict:
        met_dict = obj
    st = [number_delim.join([k,str(float(v))]) for k,v in met_dict.items()]
    return metabolite_delim.join(st)

def load_growth_medium(ex_dict, model, verbose=True):
    rxns_model = [rxn.id for rxn in model.reactions]
    for rxnid, v in ex_dict.items():
        if rxnid in rxns_model:
            model.reactions.get_by_id(rxnid).lower_bound = float(v)
        else:
            if verbose:
                print(rxnid, ', the exchange rxn, is not in the model')
    return model

def create_demand_reaction(model_in, metid):
    from copy import deepcopy
    model = deepcopy(model_in)
    rxnid = 'DM_' + metid
    model.add_reactions([cobra.Reaction(rxnid)])
    
    rxn = model.reactions.get_by_id(rxnid)
    rxn.lower_bound = 0.
    rxn.upper_bound = 0.
    rxn.reaction = metid + ' -->'
    
    return model

def test_metabolite_sink(model_in, metid):
    model = create_demand_reaction(model_in, metid)
    rxn = model.reactions.get_by_id('DM_' + metid)
    
    model.objective = {}
    rxn.objective_coefficient = 1.
    rxn.upper_bound = 1000.
    fba = model.optimize()
    
    status = True if fba.fluxes[rxn.id] > 1e-6 else False
    return status, fba

def report_mass_balance(model, chargeLim=5, verbose=True):
    '''Reports a model's mass and charge balance, except for exchange reactions and generic reactions'''
    import numpy as np
    easy_Himbal = []
    hard_Himbal = []
    check_imbal = []
    imbal_strs = dict()

    for rxn in model.reactions:
        if rxn.id[:3] == 'gen' or rxn.id[:3] == 'EX_':
            continue

        imbal = check_mass_balance_cobra(rxn, model)
        imbal_str = []
        for k,v in imbal.items():
            if type(v) in [np.int64, int] or v.is_integer():
                imbal_str.append(k + ':' + str(int(v)))
            else:
                imbal_str.append(k + ':' + str(v))
        imbal_strs[rxn.id] = ', '.join(imbal_str)

        if all([abs(i) < 1e-6 for i in imbal.values()]):
            continue

        keys = [i for i in imbal.keys() if i not in ['H', 'charge']]
        if 'H' in imbal.keys() and all([imbal[k] == 0 for k in keys]):
            if imbal['H'] == imbal['charge'] and abs(imbal['charge']) < chargeLim:
                easy_Himbal.append(rxn.id)
            elif imbal['H'] == imbal['charge'] and abs(imbal['charge']) >= chargeLim:
                hard_Himbal.append(rxn.id)
            else:
                check_imbal.append(rxn.id)
        else:
            check_imbal.append(rxn.id)
    
    if verbose:
        for rxnid in check_imbal:
            rxn = model.reactions.get_by_id(rxnid)
            if 'biomass' in rxn.id.lower() or 'normBIOM_DLTN' in rxn.id or rxn.id[:3] == 'gen':
                continue

            print(rxnid)
            print(rxn.reaction)
            for met in rxn.metabolites:
                print_list = [met.id, met.formula, str(met.charge)]
                if 'kegg.compound' in met.annotation.keys():
                    kegg = met.annotation['kegg.compound']
                else:
                    kegg = '      '
                    
                if 'formula_charge_source' in met.notes.keys():
                    fc = met.notes['formula_charge_source']
                else:
                    fc = '      '
                print('\t'.join(print_list  + [kegg, met.name, fc]))
            print(imbal_strs[rxnid])
            print()
            
    return easy_Himbal, hard_Himbal, check_imbal

def calculate_molecular_weight(formula, verbose=False):
    # Assume MW of generic group (e.g., R) to be zero 
    import sys
    sys.path.append('/home/hvdinh16/Workspace/workpy2/common/')
    from .common_params import elements_mw
    from .custom_functions import compile_elements_from_formula
    
    elem_dict = compile_elements_from_formula(formula)
    mw = 0
    for elem,coeff in elem_dict.items():
        if elem in elements_mw.keys():
            mw += coeff * elements_mw[elem]
        else:
            if verbose:
                print(elem + ' is not in the list of elements')
    
    return mw

def get_coeff_without_gam(model, biomId, gam_val):
    from collections import OrderedDict 
    import numpy as np

    atpm = ['atp_c', 'adp_c', 'h2o_c', 'pi_c', 'h_c']
    biomrxn = model.reactions.get_by_id(biomId)
    mets = atpm + sorted([met.id for met in biomrxn.metabolites.keys() if met.id not in atpm])

    bMets = OrderedDict()
    for metid in mets:
        met = model.metabolites.get_by_id(metid)
        bMets[met.id] = biomrxn.metabolites[met]

    for metid in atpm:# Exclude ATP maintainance
        bMets[metid] = bMets[metid] - np.sign(bMets[metid])*gam_val
        
    mets_rmv = []
    for k,v in bMets.items():
        if v == 0:
            mets_rmv.append(k)
            
    bMets = {k:v for k,v in bMets.items() if k not in mets_rmv}
            
    return bMets

def extract_details_from_rxnid(rxn_id):
    idsplit = rxn_id.split('-')
    tag = idsplit[0]
    rxn_base_id = idsplit[1]
    enz_id = rxn_id[(len(tag) + len(rxn_base_id) + 2):]
    rxn_dir = rxn_base_id.split('_')[-1]
    rxn_base_id = rxn_base_id[:-len(rxn_dir)-1]
    return(tag,rxn_base_id,rxn_dir,enz_id)
