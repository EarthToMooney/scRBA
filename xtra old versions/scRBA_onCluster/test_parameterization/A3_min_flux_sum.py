path_gams = '../GAMS/min_flux_sum/'
path_out = './min_flux_sum/'

import os

#### Write-out all metabolic reactions
with open('../GAMS/model/RBA_rxns.txt') as f:
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
#shutil.copy(os.path.join(path_gams, 'min_flux_violation.gms'),
#            os.path.join(path_out, 'min_flux_violation.gms'));
#shutil.copy(os.path.join(path_gams, 'soplex.opt'),
#            os.path.join(path_out, 'soplex.opt'));
            
cmds = ['cd ' + path_out,
        'module load gams',
        'gams min_flux_sum.gms  o=/dev/null']
        #'gams min_flux_sum.gms']
os.system('\n'.join(cmds))
