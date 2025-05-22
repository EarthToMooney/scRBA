* Default setup for RBA runs, including options, sets/parameters, variables, and output file names
* Reduces redundancy in the code and makes sweeping changes easier to implement
* default paths or settings specific to your model should go in phenotype_default_paths.txt or phenotype_default.txt respectively
* usage: $include <path_to_this_file>

$INLINECOM /*  */
* use default settings if not set already
$if not setGlobal gms $setGlobal gms %system.FN% 
$if not setGlobal model_root_path $setGlobal model_root_path ../ 
$if not setGlobal gms_path $setGlobal gms_path %system.FP%model/
$if not setGlobal settings_path $setGlobal settings_path %system.FP%default_paths.txt
$include %settings_path%
* default model-/phenotype-specific settings; not needed, so only loaded if it exists
$if exist %phenotype_default_paths_path% $include %phenotype_default_paths_path%
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal vmax 1e3 /* upper bound for fluxes */
$setGlobal vmin 1e-9 /* lowest positive value distinguishable from 0 by solver */
$setGlobal nscale 1e5 /* multiplier for fluxes, to reduce solver precision issues */
$setGlobal mu 0.1 /* default growth rate used as placeholder (e.g., for yields) */

* for performing bisection
** determines tolerance (i.e., stops when upper/lower bounds are negligibly far apart)
$setGlobal bi_tol 1e-5
** max number of times bisection occurs
$setGlobal max_iter 100

* Default: ≤100% of protein mass can be in mitochondria or its membrane
$setGlobal max_allowed_mito_proteome_allo_fraction 1
* Enforce part of proteome allocate to non-modeled protein
$setGlobal nonmodeled_proteome_allocation 0

* default solver settings
$include %system.FP%default_solver.gms

* remove existing modelStat file to avoid issues w/ model status detection       
file ff /%gms%.modelStat.txt/; putclose ff '';

$include %system.FP%default_sets.gms
alias (j1,j);

Parameters
S(i,j)
$include %sij_path%
dir(sm_j,j) "SM rxn, RBA rxn, and direction (-1 if RBA rxn is the reverse of SM rxn, 1 otherwise)"
$include "%sm_rxn_pairs_path%"
v_trans_UB(AA) "max molar fraction of all AAs that each AA can comprise" /#AA 1/
v_trans_LB(AA) "same as above, but min instead" /#AA 0/
NAA(j)
$include %prolen_path%
mv_exp_prot
$include %prot_mass_flux_path%
kapp(j)
$if not %ignore_kapps%==1 $include %kapp_path%
kribo(r)
$include %kribo_path%
bi_UB "upper bound for bisection" /1/
bi_LB "lower bound for bisection" /0/
max_iter /%max_iter%/
mu
prowaste_predicted(j)/
$if exist %prowaste_from_pro_with_data_path% $include %prowaste_from_pro_with_data_path%
$if exist %prowaste_from_pro_without_data_path% $include %prowaste_from_pro_without_data_path%
/
;

Variables
z, v(j), v_sum, v_trans(AA), v_trans_AA_sum;
variable table v_sm(sm_j) initial values
$include %sm_rxn_bounds_path%
;

*** FLUX LOWER AND UPPER BOUNDS (unless otherwise needed, use v_sm to alter flux bounds, to be more concise and compatible w/ other models) ***
v.lo(j) = 0; v.up(j) = %vmax%; 

* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Non-0 protein waste flux might cause error for solver by adding thousands more free variables
* Thus, protein waste should only be implemented with actual data to constrain the free variable
v.fx(j)$prowaste(j) = 0; v.up('PROWASTE-TOTALPROTEIN') = inf; v.up('PROWASTE-PROTCYT') = inf;

* Disable all biomass reactions, to avoid accidental usage
* User should activate condition-specific biomass reaction in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;

* optional model-/phenotype-specific settings 
$if exist %phenotype_default_path% $include %phenotype_default_path%
* optional task-specific phenotype file
$if exist %phenotype_path% $include %phenotype_path%
*allow bisection
mu=%mu%;

* scale all bounds by nscale, to improve solver precision; any v bounds defined after here should follow suit
v_sm.up(sm_j) = v_sm.up(sm_j) * %nscale%;
v_sm.lo(sm_j) = v_sm.lo(sm_j) * %nscale%;
v.up(j) = v.up(j) * %nscale%;
v.lo(j) = v.lo(j) * %nscale%;

Equations
Obj, Stoic, v_sum_def, v_trans_def, v_trans_AA_sum_def, v_trans_UB_fx, v_trans_LB_fx, SM_bounds, RiboCapacity, FlexibleRiboCapacity, minPro, ProwasteLim, NonModelProtAllo, MitoProtAllo, EnzCap(j,j)
;

Obj..			z =e= v('BIOSYN-PROTMODELED');
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
v_sum_def..		v_sum =e= sum(j, v(j));
v_trans_def(AA).. v_trans(AA) =e= sum(i$tRNA(AA,i),sum(j$(S(i,j) lt 0 and (sum(r,trans(r,j)) or trans_by_any_ribo(j))),v(j)*(-S(i,j))));
v_trans_AA_sum_def..	v_trans_AA_sum =e= sum(AA,v_trans(AA));
v_trans_UB_fx(AA).. v_trans(AA) =l= v_trans_AA_sum * v_trans_UB(AA);
v_trans_LB_fx(AA).. v_trans(AA) =g= v_trans_AA_sum * v_trans_LB(AA);
SM_bounds(sm_j)..	v_sm(sm_j) =e= sum(j$dir(sm_j,j), dir(sm_j,j)*v(j));
RiboCapacity(r).. sum(j$j_ribo(r,j),v(j) * kribo(r)) =g= mu * sum(j$trans(r,j), NAA(j) * v(j));
FlexibleRiboCapacity.. sum(r,sum(j$j_ribo(r,j),v(j) * kribo(r))) =g= mu * (sum(r,sum(j$trans(r,j), NAA(j) * v(j))) + sum(j$trans_by_any_ribo(j), NAA(j) * v(j)));
minPro(j)$(v.up(j) ne 0 and prowaste_predicted(j)).. v(j) * mv_exp_prot =g= v('BIOSYN-PROTTOBIO') * prowaste_predicted(j);
* Limits PROWASTE to the predicted value (if available); to allow any waste value, make the LHS >= 0. Values kept relative to total protein mass from kapp calculations for reliable scaling across different conditions.
ProwasteLim(j)$(prowaste_predicted(j)).. sum(j1$pro_syn_waste(j,j1),v(j1)) * mv_exp_prot =l= v('BIOSYN-PROTTOBIO') * prowaste_predicted(j);
NonModelProtAllo..	v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
EnzCap(j,j1)$enzload_rxn_coupling(j,j1).. v(j1) * mu =e= v(j) * kapp(j1);

*** BUILD OPTIMIZATION MODELS ***
Model 
rba /all/
fba /rba-RiboCapacity-FlexibleRiboCapacity-EnzCap/
;
rba.optfile = 1;
fba.optfile = 1;

file fi_v /%gms%.flux.txt/;
file fi_vgs /%gms%.flux_gamsscaled.txt/;
