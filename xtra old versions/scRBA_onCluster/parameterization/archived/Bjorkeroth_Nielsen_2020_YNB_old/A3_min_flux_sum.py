path_gams = '../../GAMS/parameterization/min_flux_sum/'
path_rxns_list = '../../GAMS/model/RBA_rxns.txt'
path_out = './min_flux_sum/'

import os

#### Write-out all metabolic reactions
with open(path_rxns_list) as f:
    rxns = f.read().split('\n')[1:-1]
rxns = [r[1:-1] for r in rxns]

rxns = [i for i in rxns if i.split('-')[0] == 'RXN']
rxns = ["'" + i + "'" for i in rxns]
rxns = ['/'] + rxns + ['/']

fname = os.path.join(path_out, 'rxnsMetab.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(rxns))

#### Simulation
import shutil
shutil.copy(os.path.join(path_gams, 'min_flux_sum.gms'),
            os.path.join(path_out, 'min_flux_sum.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));
            
cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_sum.gms  o=/dev/null']
        #'gams min_flux_sum.gms']
os.system('\n'.join(cmds))

#### Convert GAMS-scaled flux to actual flux
# All fluxes
fname = os.path.join(path_out, 'min_flux_sum.flux_gamsscaled.txt')
with open(fname) as f:
    fluxes = f.read().split('\n')
fluxes = [i for i in fluxes if i != '']
fluxes_new = []
for i in fluxes:
    r,vtype,val = i.split('\t')
    fluxes_new.append('\t'.join([r, vtype, str(float(val) / 1e3)]))
fname = os.path.join(path_out, 'min_flux_sum.flux.txt')
with open(fname, 'w') as f:
    f.write('\n'.join(fluxes_new))
