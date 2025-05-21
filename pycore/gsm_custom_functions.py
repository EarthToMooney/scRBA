import cobra
# based on initial file from Hoang Viet Dinh
# updated to Python 3

def test_import():
    """Test if the module imports correctly."""
    print('Import ok, common custom functions')

def secretion_rxns(model):
    """Get the list of secretions in a model (i.e., non-blocked exchanges with upper bounds ≥ 0)."""
    return [r for r in model.exchanges if r.upper_bound > 0]

def secretion_metabolites(model):
    """Get the list of secreted metabolites in a model (i.e., metabolites from rxns in 'secretion_rxns' function)."""
    return list(dict.fromkeys([m for r in secretion_rxns(model) for m in r.metabolites]))

def plot_2d_data(data, x=None,y=None,xlabel=None,ylabel=None,scale='linear',xscale=None,yscale=None,output_path=None,input_fig=None,input_ax=None,global_fig_ax=False,show_figure=True,pointsize=6,datacolor='#1f78b4',alpha=0.5,overlapping_points_in_legend=3,fontsize=12,sigfigs=2):
    """Plot 2D data.
    Defaults to scatter plot.\n
    Data is assumed to be in a pandas DataFrame with two columns (x on left) and headings as axis labels. Change as needed by updating x, y, xlabel, or ylabel.\n
    If global_fig_ax is True, the plot will be made on the global fig and ax objects. Useful for updating settings or adding additional data to the plot.\n
    overlapping_points_in_legend: number of overlapping points to show, if needed to highlight transparency. Not shown if alpha=1 (i.e., points aren't transparent).\n
    scale used by default, but setting xscale or yscale will override it for the respective axis.\n
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import linregress
	
    if x is None:
        x = data.iloc[:, 0]
    if y is None: 
        y = data.iloc[:, 1]
    if xlabel is None:
        xlabel = data.columns[0]
    if ylabel is None:
        ylabel = data.columns[1]
    scale = 'linear' if scale is None else scale
    xscale = scale if xscale is None else xscale
    yscale = scale if yscale is None else yscale

    xylim_max = max(max(x), max(y))

    if 'log' in [xscale,yscale] and (min(x) <= 0 or min(y) <= 0):
        xylim_min = min(min(x[x > 0]), min(y[y > 0]))
        print('xscale and yscale are set to ' + xscale + ' and ' + yscale + ', respectively.')
        print('Minimum value in x or y is less than or equal to 0. Setting minimum axis limit to lowest non-zero value.')
    else:
        xylim_min = min(min(x), min(y))
    
    if output_path:
        print(output_path + ':')
    # plot data
    if global_fig_ax:
        global fig, ax
    if input_fig is None or input_ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    else:
        fig = input_fig
        ax = input_ax
    # marker size is allegedly determined by how many 'M' characters can fit (picas/inch?)
    ax.scatter(x, y, color=datacolor, alpha=alpha, edgecolors=None, zorder=1, s=pointsize, lw=0)
    ax.set_xlabel(xlabel, fontsize=fontsize)  # Use xlabel as x-axis label
    ax.set_ylabel(ylabel, fontsize=fontsize)  # Use ylabel as y-axis label
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    axis_margin = 0.05
    l2 = 1.1

    for axis in ['top', 'right']:
        ax.spines[axis].set_linewidth(0)

    n = sum((x > 0) & (y > 0))
    ax.text(1, axis_margin, 'n = ' + str(n), transform=ax.transAxes, fontsize=fontsize, va='bottom', ha='right')

    ax.plot([xylim_min,xylim_max], [xylim_min,xylim_max], color='gray', ls='--', linewidth=1)

    # set minimium axis limits to lowest non-zero value if using log scale; scale up to ensure points don't get cut off by being right on the border
    ax.set_xlim([xylim_min, xylim_max*1.5])
    ax.set_ylim([xylim_min, xylim_max*1.5])

    legend_corner_left = xylim_min * axis_margin
    legend_corner_top = max(max(x), max(y))

    if alpha != 1 and overlapping_points_in_legend > 0: # no point showing overlapping points when they're opaque
        ax.text(axis_margin, 1, 'Overlapping Points', transform=ax.transAxes, fontsize=fontsize-2, va='center',ha='left')
        for i in range(overlapping_points_in_legend):
            # repeat i times
            for _ in range(0,i+1):
                ax.plot(axis_margin, 1 - axis_margin*(i + 1),transform=ax.transAxes, marker='.', color=datacolor, alpha=alpha, linestyle='None', markersize=pointsize)
            ax.text(1.5*axis_margin, 1 - axis_margin*(i + 1), i + 1, transform=ax.transAxes, fontsize=fontsize - 2, va='center', ha='left')
        # make gray box around legend
        ax.plot([axis_margin, axis_margin], [1 - axis_margin*(overlapping_points_in_legend + 1), 1], color='gray', linewidth=1)
    npx = np.array(list(x))
    npy = np.array(list(y))
    try:
        m, b, r_value, p_value, std_err = linregress(npx, npy)
        trendline = np.polyfit(npx, npy, 1)
        r_squared = r_value ** 2
        m = trendline[0]
        b = trendline[1]
        intercept_text = '' if b == 0 else ('+ ' if b > 0 else '- ') + str(abs(round(b, ndigits=sigfigs)))
        trendline_equation = f"y = {round(m, ndigits=sigfigs)}x {intercept_text}"
        print(trendline_equation)
        print(f"R^2 = {r_squared}")
        print(f"p-value = {p_value}")
    except ValueError:
        print(ValueError)

    # make x and y axes equal length
    ax.set_aspect('equal', adjustable='datalim')
    if output_path:
        plt.savefig(output_path, transparent=True, bbox_inches='tight')
        fig.savefig(output_path, transparent=True, bbox_inches='tight')
    if show_figure:
        plt.show()
    return fig, ax

def multi_plot_2d_data(datasets,rows_cols=None,x=None,y=None,xlabel=None,ylabel=None,scale='linear',xscale=None,yscale=None,output_path=None,input_fig=None,input_ax=None,global_fig_ax=False,show_figure=True,pointsize=6,datacolor='#1f78b4',alpha=0.5,overlapping_points_in_legend=3,fontsize=12,sigfigs=2):
    """Plot multiple 2D datasets on the same figure. May need some tweaking to update.\n
    rows_cols: tuple of (rows, columns) for subplots. If None, will default to (1, len(datasets)).\n"""
    import numpy as np
    import matplotlib.pyplot as plt
    letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    # note: currently doesn't support subplots with more than 26 subplots, or cases where 1 row isn't completely used
    i=0
    for d in datasets[len(datasets)-1:]: # needed to ensure multi-plot figure is created properly
        fig, ax = plot_2d_data(d,x=x,y=y,xlabel=xlabel,ylabel=ylabel,scale='log',show_figure=False,global_fig_ax=global_fig_ax,pointsize=pointsize)
    subplots_wide = rows_cols[0] if rows_cols else 1
    subplots_high = rows_cols[1] if rows_cols else len(datasets)
    figs, axes = plt.subplots(subplots_high, subplots_wide, figsize=(subplots_wide*5, subplots_high*5))
    row,col=0,0
    for d in datasets:
        # only use 2D indices if needed
        if subplots_high > 1: 
            if subplots_wide > 1:
                indices = [row,col]
            else:
                indices = row
        else:
            indices = col
        # add each plot to the figure as formatted by plot_2d_data
        fig, axes[indices] = plot_2d_data(d,x=x,y=y,xlabel=xlabel,ylabel=ylabel,scale=scale,show_figure=show_figure,input_ax=axes[indices],input_fig=fig,global_fig_ax=global_fig_ax,pointsize=pointsize)
        # add letter (a, b, c, etc.) to each subplot
        axes[indices].text(-0.05, 1.05, letters[i], transform=axes[indices].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
        i += 1
        # add subplots in reading order, looping around to next row if necessary
        if col+1 >= subplots_wide:
            row += 1
            col = 0
        else:
            col += 1
    # add title
    figs.suptitle('Comparison of SC vs RT', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('../combined-SC-vs-RT.svg')
    plt.show()

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

def remove_metabolite_from_rxn(self, metid, rxnid):
    '''Remove a metabolite from a reaction in a model.'''
    rxn = self.reactions.get_by_id(rxnid)
    met = self.metabolites.get_by_id(metid)
    if met in rxn.metabolites:
        rxn.subtract_metabolites({met: rxn.metabolites[met]})
    else:
        print('Metabolite not in reaction')
# make into a method of the Model class, allowing for easy use; relocate this to another file that's used only when cobra gets imported
# cobra.core.model.Model.remove_metabolite_from_rxn = remove_metabolite_from_rxn

def execute_command(model, model_donor, df_cmds, verbose=False):
    '''Execute commands (e.g., from a transaction log) to change a model based on a donor model.
    If commands don't require a new gene/reaction/metabolite, the model_donor can be a copy of the model.
    df_cmds: a dataframe with columns 'command', 'id', 'object_type'.'''
    attr_dict = {i:i for i in ['id', 'reaction', 'name', 'subsystem', 'lower_bound',
                               'upper_bound', 'compartment', 'notes', 'formula', 'charge']}
    attr_dict['gpr'] = 'gene_reaction_rule'
    
    comps_model = list(model.compartments.keys())
    rxns_model_donor = [rxn.id for rxn in model_donor.reactions]
    mets_model_donor = [met.id for met in model_donor.metabolites]
    # allow duplicate entries in index, so multiple commands can be executed on the same object
    df_cmds = df_cmds.reset_index(drop=True)
    
    for i in df_cmds.index:
        rxns = model.reactions
        mets = model.metabolites

        cmd = df_cmds.command[i]
        objid = df_cmds.id[i]
        obj = str(df_cmds.object_type[i])
        # if verbose:
        #     print(f'Processing command: {str(cmd)} {str(objid)} {str(obj)}')

        rxns_model = [rxn.id for rxn in model.reactions]
        mets_model = [met.id for met in model.metabolites]
        genes_model = [g.id for g in model.genes]
        
        # skip if 1st column is a comment (starts with #)
        if str(objid)[:1] == '#':
            if verbose:
                print('id of command ignored:'+objid)
            continue

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
                        if cmd1 in ['lower_bound', 'upper_bound']:
                            setattr(model.reactions.get_by_id(objid), attr_dict[cmd1], float(cmd2))
                        else:
                            setattr(model.reactions.get_by_id(objid), attr_dict[cmd1], cmd2)
                        
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
    elem_dict = OrderedDict()
    if formula in ['',None]:
        return elem_dict
    elems = re.findall('[A-Z][a-z]*', formula)
    for elem in elems:
        # handles positive and negative coefficients (negative for adjust_formula function)
        elemval = re.search('(?<=' + elem + ')[-\d.]*', formula).group()
        if elemval == '':
            elem_dict[elem] = 1.0
        else:
            elem_dict[elem] = float(elemval)
    return elem_dict

def compile_formula_from_elements(elem_dict:dict) -> str:
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
        imbal['charge'] += coeffs[met.id] * met.charge if met.charge else 0
        
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

def find_biomass_reactions(model):
    '''Finds biomass reactions in a model (SBO:0000629) and returns the properties you want (by default, ID).'''
    biomass_rxns = []
    for rxn in model.reactions:
        if 'biomass' in rxn.id.lower():
            biomass_rxns.append(rxn)
            break
        elif 'sbo' in rxn.annotation:
            if 'SBO:0000629' in rxn.annotation['sbo']:
                biomass_rxns.append(rxn)
                break
    return biomass_rxns

# check all nodes and remove any that have no path to target products
def remove_unconnected_nodes(G,nodes_to_check_connections_to):
    import networkx as nx
    continue_search = True
    while continue_search:
        nodes_removed = 0
        for node in list(G.nodes):
            if not any(nx.has_path(G,node,required_node) for required_node in nodes_to_check_connections_to):
                G.remove_node(node)
                nodes_removed += 1
        if nodes_removed == 0:
            continue_search = False
    return G

def find_pathways(model,target_products=None,flux_dict={},flux_reporting_decimals=-2,target_reactants=None,exclude_common_mets=True,paths_to_find='all_shortest_weighted',
                  node_font_size=10,edge_font_size=8,output_path=None,multigraph_mode=True):
    """Find pathways that form specified products in a model (if provided; otherwise, use all secreted metabolites) from specified reactants (if provided).

    ### Parameters
    - model : cobra.Model

    #### Optional Parameters
    - target_products : list of str, default None (converted to list of all secreted metabolites if None)
    - flux_dict : dict, default {}
        - Dict of reaction IDs as keys and fluxes as values, to help filter out irrelevant pathways. If not provided, stoichiometry is used for filtering.
        - In future versions, should ideally incorporate bounds, too.
    - flux_reporting_decimals : int, default -1
        - Number of decimal places to report for flux values in labels. If -1, doesn't report fluxes; if -2, reports as-is.
    - target_reactants : list of str, default None
        - List of metabolite IDs to use as starting points for the pathways. If not provided, checks for pathways involving all metabolites in the medium.
    - exclude_common_mets : bool, default True
        - Excludes water, lone protons, and metabolites commonly used as currency metabolites or cofactors (e.g. A(T/D/M)P, NAD(H), NADP(H), etc.), unless they're in the target products or reactants.
    - paths_to_find : str, default 'all_shortest'
        - 'all_shortest' : Find all shortest paths from target reactants to target products.
        - 'all_simple' : Find all simple paths from target reactants to target products.
    - output_path : str, default None
        - Path to save the output files to. If None, saves in current working directory.
    - multigraph_mode : bool, default True
        - If True, uses a multigraph (multiple edges allowed between the same pair of nodes in the same direction) to allow multiple edges between nodes. If False, uses a simple graph. Included only if needed for running older code from rtRBA. As of 2025-05-20, the multigraph mode makes pathway simplification less effective, so fix this.
    """
    import numpy as np, pandas as pd, networkx as nx, sys, os, random
    from networkx import Graph
    from copy import deepcopy

    medium_mets = list(i.replace('EX_','') for i in model.medium.keys() if 'EX_' in i and model.medium[i] != 0)
    # don't have specific compounds in mind, just want to see what can make the product(s)
    output_path = output_path if output_path is not None else os.getcwd() + '/'
    # check if output_path exists; if not, create it
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # force graph to use the same random seed, for easier comparisons between simulations/reruns
    seed = 1
    np.random.seed(seed)
    # Create a directed graph
    G = []
    # make new digraph
    if multigraph_mode:
        blank_graph = nx.MultiDiGraph()
    else:
        blank_graph = nx.DiGraph()
    G.append(blank_graph)

    fluxdata = pd.DataFrame(flux_dict.items(), columns=['rxn','flux'])
    rxns_used = fluxdata[fluxdata.flux != 0].rxn.values
    # find max flux for scaling
    max_flux = abs(fluxdata.flux).max()

    # only show rxns used
    # TO DO: update to use common_mets_to_exclude_from_pathway_diagrams.xlsx for common metabolites to exclude. Allow default list as well as user-provided list.
    rxns_to_consider = [r for r in model.reactions if r.id in rxns_used]
    excluded_mets = ['diphosphate','phosphate','phosphate [cytoplasm]','H2O','WATER','H','H+','PROTON','ATP','ADP','AMP','NADH','NADPH','NAD+','NADP+','NAD','NADP','NADH2','NADPH2','NADP(+)','NADP(H)','NAD(H)','NAD(P)H','NAD(P)','FAD','FADH2','FMN','FMNH2']
    # Follows same order as excluded_mets
    excluded_ids = {
        'chebi':['CHEBI:33019','CHEBI:43474','CHEBI:15377','CHEBI:24636','CHEBI:30616',"CHEBI:456216","CHEBI:456215","CHEBI:57945","CHEBI:57783","CHEBI:57540","CHEBI:58349","CHEBI:57692","CHEBI:58307","CHEBI:58210","CHEBI:57618"]
    }
    if not exclude_common_mets:
        excluded_mets = []
        # set all values in excluded_ids to empty list
        for k in excluded_ids:
            excluded_ids[k] = []

    def met_is_excluded(met):
        if met.id in excluded_mets or met.name in excluded_mets:
            return True
        elif (met.annotation.get('biocyc') is not None \
        and len(set(met.annotation.get('biocyc')).intersection(set(['META:'+i for i in excluded_mets]))) > 0) \
        or (met.annotation.get('chebi') is not None and len(set([met.annotation.get('chebi')]).intersection(set(excluded_ids['chebi']))) > 0):
            return True
        elif (met.annotation.get('chebi') is not None and len(set(met.annotation.get('chebi')).intersection(set(excluded_ids['chebi']))) > 0):
            return True
        else:
            return False

    # if target_products aren't provided, make a list of all excreted metabolites, excluding common metabolites if exclude_common_mets is True
    if target_products in [[], None, '']:
        target_product_rxns = secretion_rxns(model)
        if flux_dict != {}:
            # use flux_dict to narrow the list of secreted metabolites
            target_product_rxns = [i for i in target_product_rxns if i.id in flux_dict.keys() and flux_dict[i.id] > 0]
        target_products = list(dict.fromkeys([m.id for r in target_product_rxns for m in r.metabolites]))
        if exclude_common_mets:
            target_products = [i for i in target_products if not(met_is_excluded(model.metabolites.get_by_id(i)) or i in excluded_mets)]
        print("No target products provided; using all secreted metabolites:",target_products)
    filename_base = output_path + "pathways-making-"+str(target_products[0])+'-in-'+model.id

    if target_reactants in [[], None, '']:
        target_reactants = [i for i in medium_mets if not(met_is_excluded(model.metabolites.get_by_id(i)) or i in target_products)]
        print("No reactants provided as stopping points; using medium metabolites")
    elif type(target_reactants) == str:
        target_reactants = [target_reactants]

    for r in sorted(rxns_to_consider, key=lambda x: x.id):
        reactants,products = [r.reactants,r.products] if fluxdata[fluxdata.rxn == r.id].flux.values[0] > 0 else [r.products,r.reactants]
        for product in products:
            for reactant in reactants:
                if not(met_is_excluded(reactant)) and not(met_is_excluded(product)):
                    # uncomment line below to account for stoichiometry in the graph
                    flux = fluxdata[fluxdata.rxn == r.id].flux.values[0]*abs(r.get_coefficient(reactant))
                    label = r.id
                    if flux_reporting_decimals != -1:
                        flux_on_label = flux if flux_reporting_decimals == -2 else round(flux,flux_reporting_decimals)
                        label = r.id + '\n(' + str(flux_on_label) + ')'
                    # add fluxes as weights to edges; weight decreases with increasing abs(flux) to emphasize higher fluxes
                    # width is set to abs(flux) for visualization
                    G[0].add_edge(reactant.id,product.id,id=r.id,label=label,rounded_flux=flux_on_label,flux=flux,width=abs(flux),weight=max_flux/abs(flux))

    # set node color
    for node in G[0].nodes:
        G[0].nodes[node]['LabelGraphics'] = {'fontSize': node_font_size}
        if node in target_reactants:
            # orange
            G[0].nodes[node]['graphics'] = {'fill': '#FFA500'}
        elif node in target_products:
            # light green
            G[0].nodes[node]['graphics'] = {'fill': '#90EE90'}
        elif node.endswith('_e'):
            # make light yellow
            G[0].nodes[node]['graphics'] = {'fill': '#ffff99'}
        else:
            # make light gray
            G[0].nodes[node]['graphics'] = {'fill': '#D3D3D3'}
    # show pathways that lead to product from mets in medium
    paths = set()
    # remove all self-loops
    G.append(deepcopy(G[0]))
    G[1].remove_edges_from(nx.selfloop_edges(G[1]))
    # make new version that includes only nodes that have a path to target products
    G.append(remove_unconnected_nodes(deepcopy(G[1]),target_products))
    # make a new digraph with all nodes from G for the shortest pathways between target reactants and target products
    G.append(blank_graph)

    global mets_checked
    mets_checked = set()
    def simplest_pathway(old_graph,new_graph,target_products,possible_reactants,ptf=paths_to_find):
        global mets_checked
        for product in target_products:
            mets_checked.add(product)
            for reactant in possible_reactants:
                mets_checked.add(reactant)
                try: 
                    # find shortest path from met to target_products; include edge names in path
                    # newpaths = [nx.shortest_path(old_graph, source=reactant, target=product)]
                    if ptf == 'all_shortest_weighted':
                        newpaths = nx.all_shortest_paths(old_graph, source=reactant, target=product, weight='weight')
                    elif ptf == 'all_shortest':
                        newpaths = nx.all_shortest_paths(old_graph, source=reactant, target=product)
                    elif ptf == 'all_simple':
                        newpaths = nx.all_simple_paths(old_graph, source=reactant, target=product)
                    # find edge names in path
                    edgepath = []
                    # if any paths are found, add to set
                    # print("Paths from ",reactant," to ",product)
                    for path in newpaths:
                        for i in range(len(path)-1):
                            if old_graph.is_multigraph():
                                # Use the first available key for the edge between path[i] and path[i+1]
                                keys = list(old_graph[path[i]][path[i+1]].keys())
                                if not keys:
                                    continue  # No edge exists
                                key = keys[0]
                            else:
                                key = None  # for non-multigraphs

                            # check if new_graph already has the edge with the same label
                            already_present = False
                            if path[i] in new_graph.nodes and path[i+1] in new_graph[path[i]]:
                                if old_graph.is_multigraph():
                                    for k2 in new_graph[path[i]][path[i+1]]:
                                        # Compare labels to avoid duplicate edges
                                        if old_graph[path[i]][path[i+1]][key]['label'] == new_graph[path[i]][path[i+1]][k2]['label']:
                                            already_present = True
                                            break
                                else:
                                    if old_graph[path[i]][path[i+1]]['label'] == new_graph[path[i]][path[i+1]]['label']:
                                        already_present = True
                            if already_present:
                                continue
                            # add node attributes from old_graph if not present
                            for new_node in [path[i], path[i+1]]:
                                if new_node not in new_graph.nodes:
                                    new_graph.add_node(new_node, **old_graph.nodes[new_node])
                            # Append label to edgepath
                            if old_graph.is_multigraph():
                                edge_label = old_graph[path[i]][path[i+1]][key]['label']
                                edge_attrs = old_graph[path[i]][path[i+1]][key]
                            else:
                                edge_label = old_graph[path[i]][path[i+1]]['label']
                                edge_attrs = old_graph[path[i]][path[i+1]]
                            edgepath.append(edge_label)
                            # add edge to new_graph
                            if old_graph.is_multigraph():
                                # Let networkx assign the key automatically
                                new_edge_key = new_graph.add_edge(path[i], path[i+1])
                                # Find the last added key
                                new_keys = list(new_graph[path[i]][path[i+1]].keys())
                                new_key = new_keys[-1]
                                nx.set_edge_attributes(new_graph, {(path[i], path[i+1], new_key): edge_attrs})
                                rxn_id = edge_attrs['id']
                            else:
                                new_graph.add_edge(path[i], path[i+1])
                                nx.set_edge_attributes(new_graph, {(path[i], path[i+1]): edge_attrs})
                                rxn_id = edge_attrs['id']
                            # find all reactants in each rxn on the path
                            r = model.reactions.get_by_id(rxn_id)
                            reactants,products = [r.reactants,r.products] if fluxdata[fluxdata.rxn == r.id].flux.values[0] > 0 else [r.products,r.reactants]
                            # for reactant in reactants:
                            #     if not(met_is_excluded(reactant)) and reactant.id not in possible_reactants and reactant.id not in mets_checked:
                            #         new_graph = simplest_pathway(old_graph,new_graph,[reactant.id],possible_reactants,ptf=ptf)
                except nx.NetworkXNoPath:
                    print("No path from",reactant,"to",product)
                except nx.NodeNotFound:
                    print("Node not found:",reactant,"or",product)
                # except MaximumRecursionError:
                #     print("Maximum recursion error")
        return new_graph

    sys.setrecursionlimit(len(model.metabolites))
    G[3] = simplest_pathway(G[2],G[3],target_products,target_reactants,ptf=paths_to_find)
    G[3] = remove_unconnected_nodes(G[3],target_products)
    # export as gml
    file_path=filename_base
    for i,graph in enumerate(G):
        nx.write_gml(graph, file_path+'-G'+str(i)+".gml",stringizer=lambda x: str(x))
    for graph in [G[-1]]:
        pos=nx.spring_layout(graph)
        # nx.draw(graph, with_labels=True, node_color='yellow', edge_color='red',node_size=150)
        node_colors = [nx.get_node_attributes(graph, 'graphics')[node]['fill'] for node in graph.nodes]
        nx.draw_networkx(graph, pos, node_color=node_colors, with_labels=True, node_size=150)
        # Only draw edge labels if not a multigraph
        if not graph.is_multigraph():
            edge_labels = nx.get_edge_attributes(graph,'label')
            nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels, label_pos=0.5,font_size=edge_font_size,verticalalignment='top',rotate=False,font_color='red')
        else:
            print("Skipping edge labels: draw_networkx_edge_labels does not support multiedges.")
        # nx.draw_networkx(graph, pos=nx.spring_layout(graph), node_color=nx.get_node_attributes(graph,'graphics'), with_labels=True,node_size=150)
        # draw edge labels
        print("Nodes: ",graph.number_of_nodes())
        print("Edges: ",graph.number_of_edges())
    # nx.write_gexf(G, file_path+".gexf")
    return G

def find_pathways_without_graph(model,product_id,reactant_ids=None,flux_data=None,exclude_common_mets=True,stop_at_first=True,rxns_considered=None,pathways={}):
    """UNFINISHED: find pathways that form a specified product in a model from a specified reactant (if provided).
    \nBy default, to reduce time and improve accuracy, excludes water, lone protons, and metabolites commonly used as currency metabolites or cofactors (e.g. ATP, ADP, NADH, NADPH, etc.).
    \nIf 'flux_data' is provided, only considers reactions with non-zero fluxes.
    \nCheck for rxns that can make the product from the desired reactant, if provided. 
    \nIf none found, check all rxns that consume the reactant and produce the product, 
    \nrepeating with their respective reactants and products until the desired product is found 
    \nor all rxns are checked."""
    # pathways represented by nested dicts, 
    #   starting w/ rxns making the target product and ending w/ rxns consuming the target reactant
    target_reactants = None if reactant_ids is None else [met in model.metabolites if met.id in reactant_ids else None]
    target_prod = model.metabolites.get_by_id(product_id)
    if exclude_common_mets:
        excluded_met_names = ['H2O','WATER','H','H+','PROTON','ATP','ADP','AMP','NADH','NADPH','NAD+','NADP+','NAD','NADP','NADH2','NADPH2','NADP(+)','NADP(H)','NAD(H)','NAD(P)H','NAD(P)','FAD','FADH2','FMN','FMNH2']
        # find all mets with names matching those in the excluded list
        excluded_mets = [met for met in model.metabolites if met.name in excluded_met_names and met.id not in reactant_ids and met.id != product_id]
    # if flux_data isn't provided, use stoichiometry for filtering; otherwise, use fluxes
    if flux_data is None:
        # take only rxns that can produce target_prod (positive flux possible and target_prod in products or negative flux possible and target_prod in reactants)
        rxns_considered = model.reactions if rxns_considered is None else rxns_considered
        rxns_making_target_prod = [rxn for rxn in target_prod.reactions if rxn.upper_bound > 0 and target_prod in rxn.products or rxn.lower_bound < 0 and target_prod in rxn.reactants]
        rxns_to_check = rxns_making_target_prod
        if target_reactants is not None:
            # check first for cases where the target reactant is in the reactants of the rxn
            for rxn in rxns_making_target_prod:
                for reactant in target_reactants:
                    if rxn.upper_bound > 0 and reactant in rxn.reactants or rxn.lower_bound < 0 and reactant in rxn.products:
                        pathways.append(rxn.id)
                        # remove from rxns_to_check
                        rxns_to_check.remove(rxn)
                        if stop_at_first:
                            return pathways
                        continue
        else: # look for pathways connecting the target product to an extracellular metabolite
            for rxn in rxns_making_target_prod:
                # check if an exchange rxn makes the target product
                if rxn in model.boundary:
                    pathways.append(rxn.id)
                    # remove from rxns_to_check
                    rxns_to_check.remove(rxn)
                    if stop_at_first:
                        return pathways
                # else:


        if rxns_to_check == []:
            return pathways
        # else:
            # for rxn in rxns_to_check:

    else:
        rxns_considered = [rxn for rxn in model.reactions if flux_data[rxn.id] != 0] if rxns_considered is None else rxns_considered
        rxns_making_target_prod = [rxn for rxn in model.reactions if rxn.id in model.metabolites.get_by_id(target_prod).summary().producing_flux['reaction']]
        rxns_to_check = rxns_making_target_prod
        # if target_reactants is not None:

    
    # while a complete pathway is not found, rerun this function

    # find all products for target_prod
    for rxn in rxns_making_target_prod:
        # check if the rxn can produce target_prod
        if rxn.upper_bound > 0 and target_prod in rxn.products:
            if target_reactants in rxn.reactants:
                pathways.append(rxn.id)
                if stop_at_first:
                    break
            # else: # add the reactants of the rxn to the list of mets to check
                # for met in rxn.reactants:
        elif rxn.lower_bound < 0 and target_prod in rxn.reactants:
            if target_reactants in rxn.products:
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
                print_list = [met.id, str(met.formula), str(met.charge)]
                if 'kegg.compound' in met.annotation.keys():
                    kegg = met.annotation['kegg.compound']
                else:
                    kegg = '      '
                    
                if 'formula_charge_source' in met.notes.keys():
                    fc = met.notes['formula_charge_source']
                else:
                    fc = '      '
                print('\t'.join(print_list + [kegg, met.name, fc]))
            print(imbal_strs[rxnid])
            print()
            
    return easy_Himbal, hard_Himbal, check_imbal

def calculate_molecular_weight(formula, verbose=False):
    # Assume MW of generic group (e.g., R) to be zero 
    import sys
    from common_params import elements_mw
    
    elem_dict = compile_elements_from_formula(formula)
    mw = 0
    unknown_elems = set()
    for elem,coeff in elem_dict.items():
        if elem in elements_mw.keys():
            mw += coeff * elements_mw[elem]
        else:
            if verbose:
                unknown_elems.add(elem)
                print(elem + ' is not in the list of elements')
    
    return mw, unknown_elems

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
    """Identical to utils.py version"""
    idsplit = rxn_id.split('-')
    tag = idsplit[0]
    rxn_base_id = idsplit[1]
    enz_id = rxn_id[(len(tag) + len(rxn_base_id) + 2):]
    rxn_dir = rxn_base_id.split('_')[-1]
    rxn_base_id = rxn_base_id[:-len(rxn_dir)-1]
    return(tag,rxn_base_id,rxn_dir,enz_id)

def convert_cobra_to_gams(model, output_path='./',add_slashes=True,fwd_list_filename='RBA_rxns_FWD.txt',rev_list_filename='RBA_rxns_REV.txt',fwd_rxn_list_filename='RBA_rxns_rxnmetabolicnetworkFWD.txt',rev_rxn_list_filename='RBA_rxns_rxnmetabolicnetworkREV.txt',rxntype_filename='rxntype.txt',optstoic_rxntype_filename='rxntype_modified.txt',mets_filename='RBA_species.txt',rxns_filename='RBA_reactions.txt',sij_filename='RBA_sij.txt'):
    """Converts a COBRA model to a GAMS model, using the same format as the RBA model.
    \nUseful for working with other models/programs that don't require full RBA functionality (e.g., OptStoic)."""
    import pandas as pd
    df_eqn = pd.DataFrame(columns=['id', 'type', 'coupling_type', 'coupling_species', 'reaction'])
    c = -1
    met_list = sorted(list(set(['MET-' + i.id for i in model.metabolites if i.id != ''])))
    sij = []
    # Using Patrick's FCA code
    rxn_type_pairings = {'irrev':set(), 'reversible-fwd-half':set(), 'reversible-rev-half':set(), 'pseudoreaction':set(), 'exchange-fwd-half':set(), 'exchange-rev-half':set()} # all rxns and their respective types
    rxn_types = {'irrev': 0, 'reversible-fwd-half': 1, 'reversible-rev-half': 2, 'pseudoreaction': 3, 'exchange-fwd-half': 4, 'exchange-rev-half': 5}
    optstoic_rxn_types = {'irrev': 0, 'reversible-fwd-half': 0, 'reversible-rev-half': 0, 'pseudoreaction': 0, 'exchange-fwd-half': 4, 'exchange-rev-half': 5}
    rxn_list = []; rev_rxn_list = []; rev_list = []; fwd_rxn_list = []; fwd_list = []
    fca_list = []; optstoic_list = []

    def output_list(list):
        if add_slashes:
            return ['/'] + list + ['/']
        else:
            return list

    for rxn in model.reactions:
        # exchange rxns
        # if rxn.id[:3] == 'EX_':
        #     fwd_str = 'MET-' + [i for i in rxn.metabolites.keys()][0] + ' -->'
        #     rev_str = '--> MET-' + [i for i in rxn.metabolites.keys()][0]
        # else:
        #     fwd_str = build_stoichiometry_string(met_dict)
        met_dict = metabolites_dict_from_reaction_equation_RBA(rxn.reaction)
        met_dict = {k:v for k,v in met_dict.items() if k != ''}
        met_dict = {'MET-' + k:v for k,v in met_dict.items()}
        if rxn.upper_bound > 0:
            c += 1
            new_id = 'RXN-' + rxn.id + '_FWD-SPONT'
            # since GAMS is case-insensitive, check if new_id (regardless of case) is already in the list; add "_1", "_2", etc. if it is (check to make sure the updated version isn't in the list either)
            if new_id.lower() in [i.lower() for i in fwd_list]:
                i = 1
                while ('RXN-' + rxn.id + str(i) + '_FWD-SPONT').lower() in [i.lower() for i in fwd_list]:
                    i += 1
                print('Rxn ID already in list:', new_id, 'New ID for GAMS files:', 'RXN-' + rxn.id + str(i) + '_FWD-SPONT')
                new_id = 'RXN-' + rxn.id + str(i) + '_FWD-SPONT'
            fwd_list.append(new_id)
            if rxn.id[:3] == 'EX_':
                rxn_type_pairings['exchange-fwd-half'].add(new_id)
            else:
                rxn_type_pairings['reversible-fwd-half'].add(new_id)
            df_eqn.loc[c, 'id'] = new_id
            df_eqn.loc[c, 'type'] = 'metabolic'
            # df_eqn.loc[c, 'reaction'] = '-->' + 'MET-' + met.id
            for met, coeff in met_dict.items():
                if met == '':
                    continue
                sij.append("'"+met+"'.'"+new_id+"' "+str(coeff))
        if rxn.lower_bound < 0:
            c += 1
            met_dict = {k:-v for k,v in met_dict.items()}
            new_id = 'RXN-' + rxn.id + '_REV-SPONT'
            # since GAMS is case-insensitive, check if new_id (regardless of case) is already in the list; add "_1", "_2", etc. if it is (check to make sure the updated version isn't in the list either)
            if new_id.lower() in [i.lower() for i in rev_list]:
                i = 1
                while ('RXN-' + rxn.id + str(i) + '_REV-SPONT').lower() in [i.lower() for i in rev_list]:
                    i += 1
                print('Rxn ID already in list:', new_id, 'New ID for GAMS files:', 'RXN-' + rxn.id + str(i) + '_REV-SPONT')
                new_id = 'RXN-' + rxn.id + str(i) + '_REV-SPONT'
            rev_list.append(new_id)
            if rxn.id[:3] == 'EX_':
                rxn_type_pairings['exchange-rev-half'].add(new_id)
            else:
                rxn_type_pairings['reversible-fwd-half'].add(new_id)
            df_eqn.loc[c, 'id'] = new_id
            df_eqn.loc[c, 'type'] = 'metabolic'
            # df_eqn.loc[c, 'reaction'] = '-->' + 'MET-' + met.id
            for met, coeff in met_dict.items():
                if met == '':
                    continue
                sij.append("'"+met+"'.'"+new_id+"' "+str(coeff))
    # use rxn_type_pairings to update lists for each rxn type
    fwd_rxn_list = list(rxn_type_pairings['reversible-fwd-half']) + list(rxn_type_pairings['exchange-fwd-half'])
    fwd_rxn_list = ["'" + i + "'" for i in fwd_rxn_list if i != '/']
    with open(output_path + fwd_rxn_list_filename, 'w') as f:
        f.write('\n'.join(output_list(fwd_rxn_list)))
    rev_rxn_list = list(rxn_type_pairings['reversible-rev-half']) + list(rxn_type_pairings['exchange-rev-half'])
    rev_rxn_list = ["'" + i + "'" for i in rev_rxn_list if i != '/']
    with open(output_path + rev_rxn_list_filename, 'w') as f:
        f.write('\n'.join(output_list(rev_rxn_list)))

    fwd_list = fwd_rxn_list
    with open(output_path + fwd_list_filename, 'w') as f:
        f.write('\n'.join(output_list(fwd_list)))
    rev_list = rev_rxn_list
    with open(output_path + rev_list_filename, 'w') as f:
        f.write('\n'.join(output_list(rev_list)))

    # make FCA list from values in each key in rxn_type_pairings
    for k,v in rxn_type_pairings.items():
        fca_list += ["'" + i + "' " + str(rxn_types[k]) for i in v]
        optstoic_list += ["'" + i + "' " + str(optstoic_rxn_types[k]) for i in v]
        rxn_list += ["'" + i + "'" for i in v]
    rxn_list = sorted(list(set(rxn_list)))
    with open(output_path + rxntype_filename, 'w') as f:
        f.write('\n'.join(output_list(fca_list)))
    with open(output_path + optstoic_rxntype_filename, 'w') as f:
        f.write('\n'.join(output_list(optstoic_list)))
    with open(output_path + mets_filename, 'w') as f:
        f.write('\n'.join(output_list(met_list)))
    with open(output_path + rxns_filename, 'w') as f:
        f.write('\n'.join(output_list(rxn_list)))
    with open(output_path + sij_filename, 'w') as f:
        f.write('\n'.join(output_list(sij)))

def make_dummy_protein_stoich(length = 100, aa_standards_df = None, prot_df = None, met_dict = {}, gams_output_file = None, rxn_name='PROSYN-PROTDUMMY', dummy_metabolite_name='BIO-protdummy', mw=1):
    from collections import OrderedDict
    import pandas as pd
    if met_dict is not None:
        for met in met_dict:
            prot_st[met] = met_dict[met]
    if aa_standards_df is None:
        # find PROTEIN_amino_acid_map.txt by going up one directory until you find it in one of the directories below (check recursively)
        aa_standards_df = pd.read_csv('./build_model/input/PROTEIN_amino_acid_map.txt', sep='\t')
    prot_st = OrderedDict()
    for met in ['MET-atp_c', 'MET-h2o_c',
                'MET-adp_c', 'MET-pi_c', 'MET-h_c', 'MET-gtp_c',
                'MET-gdp_c']:
        if met not in prot_st:
            prot_st[met] = 0
    for aa in aa_standards_df.index:
        prot_st[aa_standards_df.tRNA_in[aa]] = -round(prot_df.N_AA[aa], 4)
        prot_st[aa_standards_df.tRNA_out[aa]] = round(prot_df.N_AA[aa], 4)
    for met in ['MET-atp_c', 'MET-h2o_c']:
        prot_st[met] -= 1
    for met in ['MET-adp_c', 'MET-pi_c', 'MET-h_c']:
        prot_st[met] += 1
    for met in ['MET-gtp_c', 'MET-h2o_c']:
        prot_st[met] -= 2*length
    for met in ['MET-gdp_c', 'MET-pi_c', 'MET-h_c']:
        prot_st[met] += 2*length
    prot_st[dummy_metabolite_name] = mw
    if gams_output_file:
        with open(gams_output_file, 'w') as f:
            # for each line, put the metabolite, the name, then the coefficient
            for met in prot_st:
                f.write("'" + met + "'.'" + rxn_name + "' " + str(prot_st[met]) + "\n")
    return prot_st

def read_spreadsheet(path,*args,**kwargs):
    """Make Pandas DataFrame from spreadsheet/table, automatically detecting the right file format.
    \nUnlike read_table, allows arguments found only in read_excel."""
    import pandas as pd
    if path.endswith('.xlsx'):
        return pd.read_excel(path,*args,**kwargs)
    else:
        return pd.read_table(path,*args,**kwargs)
