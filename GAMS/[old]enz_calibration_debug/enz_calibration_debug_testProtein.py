import sys
sys.path.append('../../pycore/')
from simulate import get_GAMS_modelStat

with open('./protein_list.txt') as f:
    x = f.read().split('\n')
pro_list = ['PROSYN-'+i for i in x if i != ['']]

import os
from collections import OrderedDict
prodict = OrderedDict()
for pro in pro_list:
    #print(enz)
    os.system('module load gams\n' + 'gams enz_calibration.gms --enzobj=' + pro + '  o=/dev/null')
    #os.system('module load gams\n' + 'gams enz_calibration.gms --enzobj=' + pro)
    #os.system('module load gams\n' + 'gams enz_calibration.gms --enzobj=' + pro + ' >/dev/null')
    stat = get_GAMS_modelStat('./enz_calibration.modelStat.txt')
    
    if stat == 'optimal':
        with open('./enz_calibration.enzobj.txt') as f:
            v = f.read()
        v = v.replace('\n', '')
        v = v.replace(' ', '')
        v = float(v)
        prodict[pro] = v
    else:
        prodict[pro] = 0

protext = [k+'\t'+str(v) for k,v in prodict.items()]
with open('./enz_calibration_debug.promaxflux.txt', 'w') as f:
    f.write('\n'.join(protext))
