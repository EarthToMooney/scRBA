path_gams = '../GAMS/min_flux_violation/'
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
    #rhs = "venzSyn('" + k + "') - enzfrac('" + k + "') * venzSlack"
    #rhs = "venzSyn('" + k + "') * (1 - venzDown + venzUp)"
    rhs = "venzSyn('" + k + "') * (1 - venzSlack)"
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
with open('../GAMS/model/RBA_rxns.txt') as f:
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
#shutil.copy(os.path.join(path_gams, 'min_flux_violation.gms'),
#            os.path.join(path_out, 'min_flux_violation.gms'));
#shutil.copy(os.path.join(path_gams, 'soplex.opt'),
#            os.path.join(path_out, 'soplex.opt'));
            
cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_violation.gms  o=/dev/null']
        #'gams min_flux_violation.gms']
os.system('\n'.join(cmds))
