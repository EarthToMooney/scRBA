path_gams = '../../GAMS/parameterization/min_flux_violation/'
path_rxns_list = '../../GAMS/model/RBA_rxns.txt'
path_out = './min_flux_violation/'
path_enz_level = './enz_from_proteome/enz_flux_calculation.txt'

import os

with open(path_enz_level) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']

enzsyn_dict = dict(); enzload_dict = dict()
for i in fluxes:
    enzid,v = i.split('\t')
    etype = enzid.split('-')[0]
    if etype == 'ENZSYN':
        if float(v) > 0:
            enzsyn_dict[enzid] = float(v)
    elif etype == 'ENZLOAD':
        if float(v) > 0:
            enzload_dict[enzid] = float(v)
    else:
        print('wtf')

#### Write enzyme synthesis constraints based on enzyme flux calculated from proteome
eqn_list = []; venzSyn = []

for k in enzsyn_dict.keys():
    # Constraints
    lhs = "v('" + k + "')"
    rhs = "venzSyn('" + k + "') * (1 - venzSlack) * %nscale%"
    eqn_list.append(lhs + ' =e= ' + rhs + ';')
    
    # Enzyme flux upper bounds as constants
    venzSyn.append("'" + k + "' " + str(enzsyn_dict[k]))
    
eqn_idx = ['enzSyn'+str(i) for i in range(0, len(eqn_list))]
eqn_list = ['enzSyn'+str(i)+'.. ' + eqn_list[i] for i in range(0, len(eqn_list))]

venzSyn = ['/'] + venzSyn + ['/']

fname = os.path.join(path_out, 'enzSyn_declares.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(eqn_idx))
fname = os.path.join(path_out, 'enzSyn_eqns.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(eqn_list))
fname = os.path.join(path_out, 'enzSyn_venzSyn.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(venzSyn))
    
#### List out reactions associated with enzyme with zero amount
with open(path_rxns_list) as f:
    idx = f.read().split('\n')[1:-1]
idx = [i[1:-1] for i in idx]
idx = [i for i in idx if i.split('-')[0] == 'ENZLOAD']

eqn_list = []; noEnzList = []

for i in idx:
    if i not in enzload_dict.keys():
        # Record flux variables that are in consideration
        rxnid = 'RXN-' + i.split('-', maxsplit=1)[1]
        noEnzList.append("'" + rxnid + "'")
noEnzList = ['/'] + noEnzList + ['/']

fname = os.path.join(path_out, 'noEnzRxnSet.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(noEnzList))

#### Simulation
import shutil
shutil.copy(os.path.join(path_gams, 'min_flux_violation.gms'),
            os.path.join(path_out, 'min_flux_violation.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));
            
cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_violation.gms  o=/dev/null']
        #'gams min_flux_violation.gms']
os.system('\n'.join(cmds))

#### Convert GAMS-scaled flux to actual flux
# All fluxes
fname = os.path.join(path_out, 'min_flux_violation.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
fname = os.path.join(path_out, 'min_flux_violation.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
    
# Only flux of reactions not associated with expressed enzymes
fname = os.path.join(path_out, 'min_flux_violation.flux_noEnzRxns_gamsscaled.txt')
if os.path.exists(fname):
    with open(fname) as f:
        fluxes = f.read().split('\n')
    fluxes = [i for i in fluxes if i != '']
    fluxes_new = []
    for i in fluxes:
        r,vtype,val = i.split('\t')
        fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
    fname = os.path.join(path_out, 'min_flux_violation.flux_noEnzRxns.txt')
    with open(fname, 'w') as f:
        f.write('\n'.join(fluxes_new))
