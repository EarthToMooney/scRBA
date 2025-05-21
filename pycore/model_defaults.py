# put all default model-specific parameters here (e.g., file paths). 
# Can be overridden later.

'''Naming conventions this script was designed for:
	- BiGG format generally followed unless otherwise specified
		- for stoichiometric models,
		- exchange reactions start with 'EX_'
		- compartments (e.g., for metabolites/rxns) are after the last underscore in the ID. 
			- transport rxns include the compartment being transported from (when flux is positive) followed by the compartment being transported to (when flux is negative) (e.g., CO2t_c_e transports CO2 from cytosol to extracellular space)
'''

from pycore.cobrapy_custom_extras import load_cobra_model
from cobra import Model

class RBA_options:
	"""Stores default parameters for RBA model building and simulation."""
	def __init__(self, SM_path, exchange_regex='^EX_', SM_biomass_id='BIOMASS', RBA_biomass_id='BIOMASS',
			  	 vmax=1000, compartment_regex={'default':'_[^_]+$'}):
		self.SM_path = SM_path # stoichiometric model path
		self.exchange_regex = exchange_regex
		self.SM_biomass_id = SM_biomass_id
		self.RBA_biomass_id = RBA_biomass_id
		self.vmax = vmax
		# compartment_regex is a dictionary with keys for object types (e.g., 'metabolite', 'reaction', 'enzyme')
		# for any object type not specified, the default regex will be used
		self.compartment_regex = compartment_regex
	# check if model from SM_path has multiple compartments
	def check_compartments(self):
		"""Check if the model has multiple compartments."""
		model = load_cobra_model(self.SM_path)
		self.protein_compartment = model.compartments[0]
		return len(model.compartments) > 2

# make an instance of the RBA_options class named model_defaults
RBA_defaults = RBA_options(SM_path='../build_model/input/iRhtoBatch-Rabinowitz.json')
print(RBA_defaults.SM_path)
model = RBA_defaults.load_SM()
print(model.compartments)