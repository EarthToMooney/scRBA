report_file = './binary_search_report.txt' # Text file recording binary search process
mu_tol = 1e-5; # Tolerance of upper and lower bound gap to tolerance search
maxiter = 100; # Maximum number of iteration
mu_min0 = 0; mu_max0 = 0.20; # User-set initial upper and lower bounds of mu

if mu_min0 < 0:
    print('Invalid negative value is set as the lower bound of mu. Mu is always positive. Revert to zero.')
    mu_min0 = 0

import sys
sys.path.append('../pycore/')
import json
from simulate import get_GAMS_modelStat, RBA_result

import shutil
shutil.copy('../GAMS/runRBA.gms', './runRBA.gms');
shutil.copy('../GAMS/soplex.opt', './soplex.opt');
shutil.copy('../GAMS/runRBA_GAMS_settings.txt', './runRBA_GAMS_settings.txt');

import os
# Test evaluation at zero
os.system('module load gams\n' + 'gams runRBA.gms --mu=0 o=/dev/null')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
if stat == 'infeasible':
    text = 'Model is infeasible at mu = 0, check model connectivity and other bugs. Terminate program.'
    print(text)
    with open(report_file, 'w') as f:
        f.write(text + '\n')
    quit()
    
elif stat == 'optimal':
    with open(report_file, 'w') as f:
        f.write('Pass mu = 0 test, now proceed to binary search.\n')
    
mu_min = mu_min0; mu_max = mu_max0; itercount = 0;

# Start binary search
class bcolors:
    GREEN = '\033[92m' #GREEN
    RED = '\033[91m' #RED
    RESET = '\033[0m' #RESET COLOR
    
# Evaluate min feasibility
mu = mu_min;
os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
if stat == 'infeasible':
    mu_max = mu_min
    mu_min = 0
    print(f"{bcolors.RED}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
    with open(report_file, 'a') as f:
        f.write(f"mu = {mu:.7f}, status = {stat}\n")
    
# Evaluate max infeasibility
mu = mu_max;
os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
while stat == 'optimal':
    print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}\n")
    with open(report_file, 'a') as f:
        f.write(f"mu = {mu:.7f}, status = {stat}\n")
    mu_max = 1.5*mu_max
    mu = mu_max;
    os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    
# Update min-max
mu = float(mu_min + mu_max) / 2; final_res = RBA_result();
while mu_max - mu_min > mu_tol and itercount < maxiter:
    itercount += 1
        
    os.system('module load gams\n' + 'gams runRBA.gms --mu=' + str(mu) + ' o=/dev/null')
    stat = get_GAMS_modelStat('./runRBA.modelStat.txt')
    if stat == 'optimal':
        mu_min = mu
        final_res.load_raw_flux(filepath='./runRBA.flux.txt')
        print(f"{bcolors.GREEN}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
    else:
        mu_max = mu
        print(f"{bcolors.RED}mu = {mu:.7f}, status = {stat}{bcolors.RESET}")
        with open(report_file, 'a') as f:
            f.write(f"mu = {mu:.7f}, status = {stat}\n")
        
    mu = float(mu_min + mu_max) / 2
    
# Write result
final_res.save_to_json('./RBA_result.json')
