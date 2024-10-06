import sys
sys.path.append('../../pycore/')
from simulate import get_GAMS_modelStat

with open('./enzSyn_venzSyn.txt') as f:
    x = f.read().split('\n')[1:-1]

enz_list = []
for i in x:
    enz_list.append(i.split(' ')[0][1:-1])

import os
from collections import OrderedDict
enzdict = OrderedDict()
for enz in enz_list:
    #print(enz)
    os.system('module load gams\n' + 'gams enz_calibration.gms --enzobj=' + enz + '  o=/dev/null')
    #os.system('module load gams\n' + 'gams enz_calibration.gms --enzobj=' + enz)
    #os.system('module load gams\n' + 'gams enz_calibration.gms --enzobj=' + enz + ' >/dev/null')
    stat = get_GAMS_modelStat('./enz_calibration.modelStat.txt')
    
    if stat == 'optimal':
        with open('./enz_calibration.enzobj.txt') as f:
            v = f.read()
        v = v.replace('\n', '')
        v = v.replace(' ', '')
        v = float(v)
        enzdict[enz] = v
    else:
        enzdict[enz] = 0

enztext = [k+'\t'+str(v) for k,v in enzdict.items()]
with open('./enz_calibration_debug.enzmaxflux.txt', 'w') as f:
    f.write('\n'.join(enztext))
