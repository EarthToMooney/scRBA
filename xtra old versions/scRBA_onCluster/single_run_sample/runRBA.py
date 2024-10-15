mu = 0.18 #set growth rate

import shutil
shutil.copy('../GAMS/runRBA.gms', './runRBA.gms');
shutil.copy('../GAMS/soplex.opt', './soplex.opt');
shutil.copy('../GAMS/runRBA_GAMS_settings.txt', './runRBA_GAMS_settings.txt');

import os
os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
#os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu))
