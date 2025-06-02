import pandas as pd
import cobra
from collections import OrderedDict
from copy import deepcopy
import os,sys
model_root_path = '../../'
sys.path.append(model_root_path)
from pycore.gsm_custom_functions import *

# Applies changes to models based on transaction files, to allow easier documentation of updates
# Metabolic model (COBRApy json)
modelpath = '../input/GSM_iSace1144_rba.json'
model = cobra.io.load_json_model(modelpath)

# if the input model is known, only need the transaction file path
transaction_files_and_model_versions = [
	{'file':'./transaction-file-20250530.csv'},
]
for item in transaction_files_and_model_versions:
	with model as m:
		transaction_file = item['file']
		model_output = item['model_output'] if 'model_output' in item else modelpath
		biomass_changes = item['biomass_changes'] if 'biomass_changes' in item else {}
		if 'id' in item:
			m.id = item['id']
		if 'name' in item:
			m.name = item['name']
		try:
			df_cmds = pd.read_csv(transaction_file)
			df_cmds.index = df_cmds.id.to_list()
			model_copy = deepcopy(m)
			execute_command(m, model_copy, df_cmds,verbose=True)
		except Exception as e:
			print(e)
		for attr in biomass_changes:
			setattr(m.reactions.BIOMASS, attr, biomass_changes[attr])
		cobra.io.save_json_model(m, model_output)