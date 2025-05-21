# runs kribo minimization (if user chooses to), then finds max ENZSYN fluxes for all enzymes
# update model-specific settings in kapp_options.py
from kapp_options import *

path_gams = '../../../../GAMS/parameterization/enz_from_proteome/'
path_out = './enz_from_proteome/'

run_setting_file_from = './GAMS_setting_files/enz_from_proteome_GAMS_settings.txt'
run_setting_file_to = './enz_from_proteome/enz_from_proteome_GAMS_settings.txt'

#### Create directory and copy run settings
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
shutil.copy(run_setting_file_from, run_setting_file_to);

#### Process data
with open(os.path.join(path_gams, 'pro_and_enz.txt')) as f:
    pro_list = f.read().split('\n')
pro_list = pro_list[1:-1]
pro_list = [i[1:-1] for i in pro_list]
pro_list = [i for i in pro_list if i.split('-')[0] == 'PTMPRO']

data = []; pro_data = []; pro_nodata = []
for met in pro_list:
    _,sid = met.split('-', maxsplit=1)
    if sid in df_data.index:
        pro_data.append("'PROIN-" + sid + "'")
        data.append("'PROIN-" + sid + "' " + str(df_data.loc[sid, 'vtrans (mmol/gDW/h)']))
    else:
        pro_nodata.append("'PROIN-" + sid + "'")
        
data = ['/'] + data + ['/']
pro_data = ['/'] + pro_data + ['/']
pro_nodata = ['/'] + pro_nodata + ['/']

# Write out run files
with open(os.path.join(path_out, 'proteome_data.txt'), 'w') as f:
    f.write('\n'.join(data))
with open(os.path.join(path_out, 'rxns_pro_data.txt'), 'w') as f:
    f.write('\n'.join(pro_data))
with open(os.path.join(path_out, 'rxns_pro_nodata.txt'), 'w') as f:
    f.write('\n'.join(pro_nodata))
    
#### Simulation
shutil.copy(os.path.join(path_gams, 'enz_from_proteome.gms'),
            os.path.join(path_out, 'enz_from_proteome.gms'));
shutil.copy(os.path.join(path_gams, 'soplex.opt'),
            os.path.join(path_out, 'soplex.opt'));

# if find_kribo:

cmds = ['cd ' + path_out,
        'module load gams',
        'gams enz_from_proteome.gms' + output_redirect_str]
os.system('\n'.join(cmds))
fname = os.path.join(path_out, 'enz_from_proteome.modelStat.txt')
stop_if_run_failed(fname)
