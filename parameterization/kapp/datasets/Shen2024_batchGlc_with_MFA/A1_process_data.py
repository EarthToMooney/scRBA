# update model-specific settings in kapp_options.py
from kapp_options import *
from copy import deepcopy

def convert_protein_abundance(abundance, output_units):
    if proteomics_units == output_units:
        return abundance
    elif proteomics_units == 'g/g protein' and output_units == 'g/gDW':
        return abundance * ptot
    elif proteomics_units == 'g/gDW' and output_units == 'g/g protein':
        return abundance / ptot
    else:
        raise ValueError("Invalid units specified. Use 'g/g protein' or 'g/gDW'.")

# load SM model
model = load_cobra_model(sm_path)

aa_map = pd.read_csv(aa_mapping_path, sep='\t')
aa_dict = dict(zip(aa_map['aa_abbv'], aa_map['MW']))

# Mass fraction of proteome that's not explicitly modeled. Set to 0 if recalculate_nonmodeled_proteome_allocation = True.
nonmodel_protein_data_path = './nonmodeled_proteins.json'
# if path exists
if os.path.exists(nonmodel_protein_data_path): 
    with open(nonmodel_protein_data_path, 'r') as f:
        # load from file
        nonmodel_proteins = json.load(f)
else:
    nonmodel_proteins = []
recalculate_nonmodeled_proteome_allocation = False
# Max mass fraction of modeled proteome that's inside the mitochondria. Set to 0 by default and automatically adjusted if recalculate_mito_proteome_allocation = True.
max_allowed_mito_proteome_allo_fraction = 0
recalculate_mito_proteome_allocation = True
ATP_cost_of_translation = 0 # mmol ATP/(gDW*h); calculated from data if 0

search_uniprot_for_nonmodeled_sequences = False
uniprot_url = 'https://rest.uniprot.org/uniprotkb/'

# flux data (e.g., from MFA) is optional
# runs only if flux data file exists
try:
    df_flux = read_spreadsheet(flux_path)
    if 'fluxes_to_ignore' in locals():
        df_flux = df_flux[~df_flux['id'].isin(fluxes_to_ignore)]
    # make a copy of the model to use for data processing (for comparing rxns)
    model_from_data = deepcopy(model)
    model_from_data.remove_reactions(model_from_data.reactions)
except FileNotFoundError:
    # make empty dataframe
    df_flux = pd.DataFrame(columns=['id', col_LB, col_UB])
df_flux.index = df_flux['id'].to_list()
# make GAMS file for flux data

# load path_model + 'SM_rxn_bounds.txt' as df (text between spaces represents the reaction ID, lower bound, and upper bound)
df_flux_bounds = pd.read_csv(path_model + 'SM_rxn_bounds.txt', sep=r'\s+')
df_flux_bounds.index.name = 'id'
rxns_to_review = dict() # key: reaction ID, value: list of comments to review

# whether to compare stoichiometry of flux data with model stoichiometry
ignore_compartments_in_flux_data_stoich = True # if True, ignores compartments in flux data stoichiometries when comparing with model stoichiometries
if 'reaction' in df_flux.columns:
    check_stoich = True
else:
    print("No 'reaction' column in flux data; skipping stoichiometry check.")
    check_stoich = False

# write ID, then lower bound (-vmax if none) and upper bound (vmax if none)
with open('./v_exp_lb.txt', 'w') as f, open('./v_exp_ub.txt', 'w') as f2, open('./sm_j_lumped.txt', 'w') as f3, open('./sm_j_lumped_mappings.txt', 'w') as f4:
    f.write('/\n'); f2.write('/\n'); f3.write("/\n'placeholder'\n"); f4.write('/\n*Uncomment lines to use them (e.g., if v=0 for \'DGAT\', which represents \'DGAT_c\' and \'DGAT_m\', line should be "\'DGAT\'.(\'DGAT_c\',\'DGAT_m\') 0")\n')
    for i in df_flux.index:
        id = "'"+i+"'"
        # check if rxn in GSM_rxn_ids.txt
        if id not in open(gsm_rxn_ids_path).read():
            print(f"'{i}' not in {gsm_rxn_ids_path}; adding to sm_j_lumped_mappings.txt")
            f3.write(f"'{i}'\n")
            f4.write(f"*'{i}'.('','') {df_flux.loc[i, col_LB]}\n") # if needed in the future, add support for upper bounds too
            continue
        else:
            # if flux outside of df_flux_bounds, print warning
            if id in df_flux_bounds.index:
                if df_flux.loc[i, col_LB] < df_flux_bounds.loc[id, 'lo'] or df_flux.loc[i, col_UB] > df_flux_bounds.loc[id, 'up']:
                    print(f"Warning: Flux bounds for '{i}' ({df_flux.loc[i, col_LB]}, {df_flux.loc[i, col_UB]}) are outside of the model bounds ({df_flux_bounds.loc[id, 'lo']}, {df_flux_bounds.loc[id, 'up']}) in {path_model + 'SM_rxn_bounds.txt'}.")
                elif check_stoich: # ensure that direction/stoichiometry don't conflict with model stoichiometry
                    rxn = cobra.Reaction(i)
                    if i not in model_from_data.reactions:
                        model_from_data.add_reactions([rxn])
                    model_from_data.reactions.get_by_id(i).reaction = df_flux.loc[i, 'reaction']
                    stoich_in_data_raw = dict(sorted({k.id:v for k,v in model_from_data.reactions.get_by_id(i).metabolites.items()}.items()))
                    stoich_in_model_raw = dict(sorted({k.id:v for k,v in model.reactions.get_by_id(i).metabolites.items()}.items()))

                    if ignore_compartments_in_flux_data_stoich:
                        # remove compartment info from flux data stoichiometry
                        stoich_in_data = {k.id.rsplit('_',1)[0]:v for k,v in model_from_data.reactions.get_by_id(i).metabolites.items()}
                        stoich_in_model = {k.id.rsplit('_',1)[0]:v for k,v in model.reactions.get_by_id(i).metabolites.items()}
                    else:
                        stoich_in_data = stoich_in_data_raw
                        stoich_in_model = stoich_in_model_raw
                    stoich_in_data = dict(sorted(stoich_in_data.items()))
                    stoich_in_model = dict(sorted(stoich_in_model.items()))
                    if stoich_in_model != stoich_in_data and {k.id:v for k,v in model.reactions.get_by_id(i).metabolites.items()} != {k.id:v for k,v in model_from_data.reactions.get_by_id(i).metabolites.items()}:
                        # if the only difference is the direction (i.e., signs of the coefficients), then make the signs match
                        if {k:-v for k,v in stoich_in_model.items()} == stoich_in_data:
                            err = 'Model and data have different directions'
                            df_flux.loc[i, 'reaction'] = model.reactions.get_by_id(i).reaction
                            # invert the signs in the flux data, and switch the upper and lower bounds
                            df_flux.loc[i, col_LB], df_flux.loc[i, col_UB] = -df_flux.loc[i, col_UB], -df_flux.loc[i, col_LB]
                        else:
                            err = 'Model and data have different stoichiometries'
                            # print(f"Warning: Stoichiometry for '{i}' in flux data does not match the model stoichiometry: {stoich_in_data} vs {stoich_in_model}.")
                        rxns_to_review.setdefault(err, dict()).setdefault(i, []).append(f"(model: {model.reactions.get_by_id(i).reaction}, flux data: {df_flux.loc[i, 'reaction']})")
            else:
                print(f"Warning: Flux bounds for '{i}' not found in {path_model + 'SM_rxn_bounds.txt'}.")
            if pd.isnull(df_flux.loc[i, col_LB]):
                lb = -vmax
            else:
                lb = df_flux.loc[i, col_LB]
            if pd.isnull(df_flux.loc[i, col_UB]):
                ub = vmax
            else:
                ub = df_flux.loc[i, col_UB]
            f.write("'"+i+"'" + ' ' + str(lb) + '\n')
            f2.write("'"+i+"'" + ' ' + str(ub) + '\n')
    if rxns_to_review:
        for err, rxns in rxns_to_review.items():
            print(f"{err}:")
            for rxn_id, details in rxns.items():
                print(f"  {rxn_id}: {', '.join(details)}")

    # save to fluxomics.tsv if df_flux is not empty
    if not df_flux.empty:
        df_flux.to_csv('./fluxomics.tsv', sep='\t', index=True, header=True)
    f.write('/'); f2.write('/'); f3.write('/\n'); f4.write('/\n')


# Load protein
df_prot_raw = read_spreadsheet(prot_path)
df_prot = df_prot_raw.copy()
df_prot.index = df_prot.id.to_list()
# Strip compartment
df_prot.index = [i.split('_')[0] if '_' in i else i for i in df_prot.index]
df_prot['id'] = df_prot.index.to_list()
df_prot = df_prot.drop_duplicates(subset=['id'])
# Protein copy selector: start with empty file
# per Hoang: "protein_copies_selector.txt" is a manually written file. It should be empty for you in the first run. Then, you have to check your kapp calculation to see if you need to add any entries. For example, an enzyme ILV2-ILV6 complex and ILV2 are functional, but I choose to assign estimated kapp value to ILV2-6 complex only (which I think is the primary complex for catalysis, otherwise it doesn't make sense to me why ILV2-6 complex is even needed). There will be many cases like this where you have to write entries to the "protein_copies_selector.txt" file
try:
    df_select = pd.read_csv('./input/protein_copies_selector.txt', sep='\t')
    df_select.index = df_select.gene_src.to_list()
except: 
    df_select = pd.DataFrame(columns=['gene_src', 'selected_compartmental_copy'])

# Find all ribosomal proteins
with open(path_gams + 'model/ribo_subunits.json', 'r') as f:
    ribo_dict = json.load(f)

#### HANDLE MISSING MEASUREMENTS FOR SUBUNIT COMPONENT OF HETEROMERIC ENZYMES
# E.g., missing subunit measurements for ATP synthase complex
# Stoichiometry
df_eqn = read_spreadsheet(model_xlsx_path)
df_eqn.index = df_eqn.id.to_list()

headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
}

for aa in aa_dict:
    dummy_protein['AA abundances'][aa] = 0
total_dummy_abundance_per_mw = 0
# find nonmodeled proteome allocation
if recalculate_nonmodeled_proteome_allocation:
    for i in df_raw.index:
        if i not in df_prot.index:
            seq = ''
            conc = 0
            mw = 0
            # print(i, "not in df_prot")
            # add to dummy protein and nonmodeled proteome allocation calculations
            c_avg = np.mean(df_raw.loc[i, cols_data])
            avg_abundance = convert_protein_abundance(c_avg, 'g/g protein')
            nonmodeled_proteome_allocation += avg_abundance
            # check nonmodeled_proteins for sequence, MW, and conc
            for p in nonmodel_proteins:
                if p['id'] == i:
                    seq = p['sequence'].replace('*','')
                    conc = p['conc (g/gDW)']
                    mw = p['MW (g/mmol)']
                    break
            # consider scraping uniprot for sequence or using API if they have one
            if search_uniprot_for_nonmodeled_sequences and not seq:
                print('searching uniprot for', df_raw.loc[i, uniprot_col])
                # search uniprot for protein sequence
                url = uniprot_url + df_raw.loc[i, uniprot_col] + '?format=json'
                # get response, convert to dict
                # response = requests.get(url).json()
                response = requests.get(url,headers=headers).json()
                # print('response:',response)
                # find "sequence" key
                if 'sequence' in response:
                    # add sequence to dummy protein
                    seq = response['sequence']['value'].replace('*','')
                    gene_src = ''
                    for gene in response['genes']:
                        # print('gene:',gene)
                        if 'geneName' in gene.keys():
                            gene_src = gene['geneName']['value']
                            # print('gene_src:',gene_src)
                            break
                    # convert this excel formula into a method of determining protein mass: =SUMPRODUCT((LEN([@sequence])-LEN(SUBSTITUTE([@sequence],{"A";"C";"D";"E";"F";"G";"H";"I";"K";"L";"M";"N";"P";"Q";"R";"S";"T";"V";"W";"Y"},""))),{72.08;104.14;115.08;129.11;148.17;58.05;138.14;114.16;130.18;114.16;132.2;115.1;98.12;129.13;158.19;88.08;102.1;100.13;187.21;164.17})/1000
                    mw = sum([len(seq) - len(seq.replace(aa, '')) for aa in aa_dict] * np.array(list(aa_dict.values()))) / 1000
                    # abundance / wt.
                    conc = convert_protein_abundance(c_avg, 'g/gDW')
                    # with open('./nonmodeled_proteins.json', 'a') as f:
                    #     f.write(str({'id':i,'URL':url,'sequence':seq,'MW (g/mmol)':mw,'conc (g/gDW)':conc}))
                    nonmodel_proteins.append({'id':i,'gene_src':gene_src,'URL':url,'sequence':seq,'MW (g/mmol)':mw,'conc (g/gDW)':conc})
            if seq:
                total_dummy_abundance_per_mw += conc / mw
                for aa in aa_dict:
                    # find fraction of amino acid in sequence, multiply by abundance / MW of protein
                    dummy_protein['AA abundances'][aa] += (((len(seq) - len(seq.replace(aa, ''))) / len(seq)) * conc / mw)
        else:
            seq = df_prot.loc[i, 'sequence']
            conc = convert_protein_abundance(c_avg, 'g/gDW')
            mw = df_prot.loc[i, 'MW (g/mmol)']
        if seq == '' or conc == 0 or mw == 0:
            print(i, conc, mw, seq)
        if conc and mw and seq != '':
            # print(i, conc, mw, seq.replace('*',''))
            ATP_cost_of_translation += (conc * ptot * ((len(seq.replace('*','')) * 2) + 1) / mw)
    # save nonmodeled protein info to JSON
    with open(nonmodel_protein_data_path, 'w') as f:
        json.dump(nonmodel_proteins, f)
    for p in nonmodel_proteins:
        # check if p['id'] is in df_prot.index
        if p['id'] not in df_prot_raw.index:
            # add to bottom of prot_path file, using "id" value in both "id" and "uniprot" columns
            df_prot_raw = pd.concat([df_prot_raw, pd.DataFrame({'id': [p['id']], 'gene_src': [p['gene_src']], 'name': [p['id']], 'uniprot': [p['id']], 'subloc_assigned': ['unknown'], 'cofactor_comments': ['Unknown; protein added automatically to help with fitting translation data'], 'MW (g/mmol)': [p['MW (g/mmol)']], 'sequence': [p['sequence']], 'status': ['forProteomicsOnly'], 'translation_loc': ['unknown']})], ignore_index=True)
            df_prot = pd.concat([df_prot, pd.DataFrame({'id': [p['id']], 'gene_src': [p['gene_src']], 'name': [p['id']], 'uniprot': [p['id']], 'subloc_assigned': ['unknown'], 'cofactor_comments': ['Unknown; protein added automatically to help with fitting translation data'], 'MW (g/mmol)': [p['MW (g/mmol)']], 'sequence': [p['sequence']], 'status': ['forProteomicsOnly'], 'translation_loc': ['unknown']})], ignore_index=True)
        else:
            print(p['id'], "already in df_prot")
    # save file
    if prot_path.endswith('.xlsx'):
        df_prot_raw.to_excel(unknown_prot_path, index=None)
    else:
        df_prot_raw.to_csv(unknown_prot_path, index=None)
    # max_allowed_mito_proteome_allo_fraction = 1 - nonmodeled_proteome_allocation

    # find median length of nonmodeled proteins
    dummy_protein['length'] = np.median([len(p['sequence']) for p in nonmodel_proteins])
    # divide all amino acid abundances by total_dummy_abundance_per_mw and multiply by length
    for aa in dummy_protein['AA abundances']:
        if total_dummy_abundance_per_mw != 0:
            dummy_protein['AA abundances'][aa] /= total_dummy_abundance_per_mw
        dummy_protein['AA abundances'][aa] *= dummy_protein['length']

    dummy_protein['MW (g/mmol)'] = sum([dummy_protein['AA abundances'][aa] * aa_dict[aa] for aa in aa_dict]) / 1000
    # convert into dataframe with 'AA abundances' keys as aa_abbv column, and values as N_AA column
    dummy_protein_df = pd.DataFrame(columns=['aa_abbv', 'N_AA'])
    # aa_abbv column from 'AA abundances' keys
    dummy_protein_df['aa_abbv'] = dummy_protein['AA abundances'].keys()
    # N_AA column from 'AA abundances' values
    dummy_protein_df['N_AA'] = dummy_protein['AA abundances'].values()
    # print(dummy_protein)
    # make rxn equation for dummy protein
    print('dummy protein:',dummy_protein)
    make_dummy_protein_stoich(prot_df=dummy_protein_df, length=dummy_protein['length'], gams_output_file='./RBA_sij_for_kapps.txt', aa_standards_df=aa_map, mw=dummy_protein['MW (g/mmol)'], rxn_name=dummy_protein['id'])
    # copy over all lines (except ones with dummy_protein['id'] in them) from sij_path to RBA_sij_for_kapps.txt
    lines = []
    with open(sij_path, 'r') as f, open('./RBA_sij_for_kapps.txt', 'r') as f2:
        for line in f2:
            if line not in ['/','/\n']:
                lines.append(line)
        for line in f:
            if ".'"+dummy_protein['id']+"'" not in line and line not in ['/','/\n']:
                lines.append(line)
    with open('./RBA_sij_for_kapps.txt', 'w') as f:
        # add slashes as first and last lines
        f.write('/\n')
        for line in lines:
            f.write(line)
        f.write('/')

# Process data
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)', 'type', 'conc (g/gDW)', 'vtrans (mmol/gDW/h)']
idx = [i for i in df_prot.index if i in df_raw.index]

df_data = pd.DataFrame(columns=cols, index=idx)
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)']
df_data.loc[idx, cols] = df_prot.loc[idx, cols]

protein_categories = dict() # key: protein id, value: set of categories to track (e.g., mitochondrial, modeled)
for i in df_data.index:
    protein_categories[i] = set()
    data = df_raw.loc[i, cols_data].values.tolist()
    # remove rows containing the column name; print out which were removed
    # for d in data:
    #     if type(d) == str:
    #         print(i, 'dropped from data due to being a string:', d)
    #         data = data.drop(d)
    # for d in data:
    #     if type(d) != ('int' or 'float'): 
    
    # Filter out null values
    data = [c for c in data if pd.isnull(c) == False]
    # print(data)
    #if data == []:
    #    df_data.loc[i, 'conc (g/gDW)'] = 0
    #    df_data.loc[i, 'vtrans (mmol/gDW/h)'] = 0
    if data != []:
        # print(data)
        c_list = []
        for c in data:
            if type(c) != 'int' and type(c) != 'float': 
                continue
            else:
                c_list.append(c)
        # Average out the data for each protein (accounts for replicate experiments, if done)
        c_avg = np.mean([data])
        # c_avg = np.mean([c_list])
        # c_avg = 
        # # print(c_avg)
        # # c_avg = data
        mw = df_prot.loc[i, 'MW (g/mmol)']
        df_data.loc[i, 'c_avg'] = c_avg
        df_data.loc[i, 'conc (g/gDW)'] = convert_protein_abundance(c_avg, 'g/gDW')
        df_data.loc[i, 'vtrans (mmol/gDW/h)'] = mu * df_data.loc[i, 'conc (g/gDW)'] / mw
        df_data.loc[i, 'type'] = 'truedata_enz'
        
        if i in ribo_dict.keys():
            df_data.loc[i, 'type'] = 'truedata_ribo'
            if 'm' in ribo_dict[i]:
                protein_categories[i].add('can be in mitochondria')

# Reindex - incorporate info from protein copy selector
idx = [df_select.selected_compartmental_copy[i] if i in df_select.index        else i for i in df_data.index]
df_data.index = idx
df_data['id'] = df_data.index.to_list()

# Clean out NaN rows
df_data = df_data[df_data['conc (g/gDW)'].isnull() == False]

# Gap-fill data
# Load protein
df_prot = read_spreadsheet(prot_path)
df_prot.index = df_prot.id.to_list()

# for each row, check if it's in selected_compartmental_copy
# if it is, replace the index with the selected_compartmental_copy value
# if not, then add all compartment-specific copies of it to the output file

# copy df_data
iter = 0
sumlimits_proin = []
sumlimits = []
sumlimits_pro_set = []
pro_data = []
df_data_copy = df_data.copy()
for i in df_data.index:
    if i in df_select.index:
        df_data_copy.loc[i, 'id'] = df_select.loc[i, 'selected_compartmental_copy']
        print('selected compartmental copy:',i, df_select.loc[i, 'selected_compartmental_copy'])
    else: 
        if i.split('_')[0] in ['m','mm']:
            protein_categories[i].add('can be in mitochondria')
        # comment out lines below to equally distribute protein abundance among all matches
        conc = df_data.loc[i, 'conc (g/gDW)']
        vtrans = df_data.loc[i, 'vtrans (mmol/gDW/h)']
        if i in df_prot['gene_src'].values:
            matches = list(df_prot.loc[df_prot['gene_src'] == i].iterrows())
            #conc = df_data.loc[i, 'conc (g/gDW)'] / len(matches)
            #vtrans = df_data.loc[i, 'vtrans (mmol/gDW/h)'] / len(matches)
            # make list for all names of protein copies
            allcopies = []
            for index, row in df_prot.loc[df_prot['gene_src'] == i].iterrows():
                if row['gene_src'] == i:
                    new_row = df_data.loc[i].copy()
                    new_row['id'] = row['id']
                    if row['subloc_assigned'] in ['mm','m']: # if the protein is potentially mitochondrial
                        protein_categories[i].add('can be in mitochondria')
                    allcopies.append(row['id'])
                    # divide conc and vtrans by the number of matches, evenly splitting the protein abundance
                    # new_row['conc (g/gDW)'] = new_row['conc (g/gDW)'] / len(matches)
                    # new_row['vtrans (mmol/gDW/h)'] = new_row['vtrans (mmol/gDW/h)'] / len(matches)
                    new_row['conc (g/gDW)'] = conc
                    new_row['vtrans (mmol/gDW/h)'] = vtrans
                    df_data_copy.loc[len(df_data_copy)] = new_row
            iter += 1
            # new_row = df_data.loc[i]
            # new_row['id'] = df_prot.loc[df_prot['gene_src'] == i, 'id'].values[0]
            # df_data_copy = df_data_copy.concat(new_row, ignore_index=True)
            # for index, row in df_prot.iterrows():
            #     if row['gene_src'] == i:
            #         new_row = df_data_copy.loc[i]
            #         new_row['id'] = row['id']
            #         df_data_copy = df_data_copy.append(new_row, ignore_index=True)
        else:
            iter += 1
            allcopies = [i]
        sumlimits_pro_set.append("'" + i + "'")
        pro_data.append(f"'{i}' {vtrans}")
        if vtrans == 0 and allow_trans_when_measurement_is_0:
            sumlimits_proin.append("Equation prosum" + str(iter) + "; prosum" + str(iter) + ".. " + " + ".join(["v('PROIN-" + copy + "')" for copy in allcopies]) + " =l= " + str(vtrans) + "*1e6;")
            sumlimits.append("Equation prosum" + str(iter) + "; prosum" + str(iter) + ".. " + " + ".join(["v('PROSYN-" + copy + "')" for copy in allcopies]) + " =e= prosynSlackUB('" + i + "');")
        else:
            sumlimits_proin.append("Equation prosum" + str(iter) + "; prosum" + str(iter) + ".. " + " + ".join(["v('PROIN-" + copy + "')" for copy in allcopies]) + " =l= " + str(vtrans) + "*1e6;")
            sumlimits.append("Equation prosum" + str(iter) + "; prosum" + str(iter) + ".. " + " + ".join(["v('PROSYN-" + copy + "')" for copy in allcopies]) + " =e= " + str(vtrans) + " * %nscale% * (1 - prosynSlackLB('" + i + "') + prosynSlackUB('" + i + "'));")
    if recalculate_mito_proteome_allocation and 'can be in mitochondria' in protein_categories[i]:
        max_allowed_mito_proteome_allo_fraction += df_data.loc[i,'c_avg']
# remove all duplicate rows
df_data_copy = df_data_copy.drop_duplicates(subset=['id'], keep='first').sort_values('id')
df_data_copy_filtered = df_data_copy.copy()

# create constraints for protein abundance
with open('./unique_proteins_no_locations.txt', 'w') as f:
    f.write("\n".join(['/'] + sumlimits_pro_set + ['/']))
with open('./protein_abundance_constraints.txt', 'w') as f:
    for limit in sumlimits_proin:
        f.write(limit + '\n')
with open('./prosyn_abundance_constraints.txt', 'w') as f:
    for limit in sumlimits:
        f.write(limit + '\n')
with open('./v_exp_pro.txt','w') as f:
    f.write('\n'.join(['/'] + sorted(pro_data) + ['/']))

# if any row has an "id" value not in the "id" values of df_prot, then print that row
errors = []
for i in df_data_copy.index:
    if df_data_copy.loc[i, 'id'] not in df_prot.index:
        # important for removing redundant rows (e.g., ones w/o locations for proteins that always have them)
        errors.append(str(i))
        # remove row
        df_data_copy_filtered = df_data_copy_filtered.drop(i)
df_data_copy_filtered.index = df_data_copy_filtered['id'].to_list()
df_data_copy_filtered

idx_enzsyn = df_eqn[df_eqn.id.str.contains('ENZSYN-')].index
cols = ['id', 'name', 'uniprot', 'MW (g/mmol)']
# print(df_prot)
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
    # print(met_dict)
    in_data = set(met_dict) & set(df_data_copy_filtered.index)
    if len(in_data) > 0.5 and len(in_data) < len(met_dict): # i.e., some but not all subunits are measured
        #print(i, len(in_data), len(met_dict), ','.join(in_data))
        vmin = min([df_data_copy_filtered.loc[k, 'vtrans (mmol/gDW/h)'] / met_dict[k] for k in met_dict.keys() if k in in_data])
        cmin = min([df_data_copy_filtered.loc[k, 'conc (g/gDW)'] / met_dict[k] for k in met_dict.keys() if k in in_data])
        for k in met_dict.keys():
            # print(df_data)
            if k not in in_data:
                # print("k=",k)
                df_data_copy_filtered.loc[k, cols] = df_prot.loc[k, cols]
                df_data_copy_filtered.loc[k, 'conc (g/gDW)'] = cmin * met_dict[k]
                df_data_copy_filtered.loc[k, 'vtrans (mmol/gDW/h)'] = vmin * met_dict[k]
                df_data_copy_filtered.loc[k, 'type'] = 'gapfill_subunit'

df_data_copy_filtered.to_csv(path_data, index=None, sep='\t')

# print(df_data_copy_filtered[df_data_copy_filtered.duplicated(subset='uniprot', keep=False)].sort_values('uniprot'))

# show all gapfill_subunit rows
if df_data_copy_filtered[df_data_copy_filtered['id'] == 'gapfill_subunit'].shape[0] > 0:
    print('gapfill_subunit rows:')
    print(df_data_copy_filtered[df_data_copy_filtered['id'] == 'gapfill_subunit'])

# if any row has an "id" value not in the "id" values of df_prot, then print that row
if errors:
    error_message = "Protein IDs not in "+prot_path+":\n" + ", ".join(errors)
    # raise ValueError(error_message)
    print(error_message)

# print nonmodeled proteome allocation
print('nonmodeled proteome allocation:', nonmodeled_proteome_allocation)
# print max allowed mito proteome allocation
print('max allowed mito proteome allocation:', max_allowed_mito_proteome_allo_fraction)
print('ATP costs of translation:',ATP_cost_of_translation)
