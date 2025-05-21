sets
i
$include %species_path%
pro
$include %pro_path%
AA
$include %aa_path%
tRNA(AA,i)
$include %tRNA_path%
j
$include %rxns_path%
g
$include %genes_path%
expression(g,j)
$include %expression_path%
sm_j "stoichiometric model (SM) reactions"
$include %sm_rxns_path%
prosyn(j) 
$include %prosyn_path%
pro_prosyn(pro,j)
$include %proteins_and_locations_path%
pro_syn_waste(j,j)
$include %model_root_path%build_model/model/pro_syn_waste.txt
r 
$include %ribosomes_path%
pro_ribo(r,pro)
$include %gms_path%pro_ribo_subunits.txt
j_ribo(r,j)
$include %j_ribo_path%
trans(r,j)
$include %translation_path%
trans_by_any_ribo(j)
$include %trans_by_any_ribo_path%
prowaste(j)
$include %prowaste_path%
uptake(j) "uptake rxns"
$include %uptake_path%
media(j) "allowed uptake rxns based on simulated media conditions"
$include %media_path%
rxns_biomass(j)
$include %biomass_path%
enzload_rxn_coupling(j,j)
$include "%enzload_rxn_coupling_path%"
;
