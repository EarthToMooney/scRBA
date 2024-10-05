*** Minimize protein production flux while adjusting kapps ***
*       Authors: Eric Mooney
***********************************************************************************

$INLINECOM /*  */
$include "../GAMS_setting_files/test_kapp_GAMS_settings.txt"
$setGlobal nscale 1e5
* max fluxes allowed, and min fluxes deemed significant enough to report
$setGlobal vmax 1e3
$setGlobal vmin 0 

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1

* Default value for ribosome efficiency
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

$setGlobal nonmodeled_proteome_allocation 0
$setGlobal venzSlackAllow 0
$setGlobal prosynSlackAllow 0
* small value needed to ensure sequential problems aren't infeasible due to rounding errors
$setGlobal epsilon 1e-7

options
	LP = cplex /*Solver selection*/
	limrow = 0 /*number of equations listed, 0 is suppresed*/
	limcol = 0 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 10000 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
pro
$include "%unique_protein_set_path%"
prowaste(j)
$include "%prowaste_path%"
gsm_j /* list of GSM model rxns */
$include "%gsm_rxns_path%"
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
unknown_ribo_translation(j)
$include "%unknown_ribo_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_with_no_prodata(j) /*enzymatic rxns without proteomics data for all their enz subunits*/
$include "%rxns_with_no_prodata_path%"
prodata_set(j)
$include "%proteome_data_set_path%"
rxns_metab(j)
$include "%rxns_metab_path%"
rxns_biomass(j) /*list of biomass rxns*/
$include "%biomass_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
pro_val(j)
$include "%proteome_data_path%"
dir(gsm_j,j) /* lists GSM rxn, RBA rxn, and direction (-1 if RBA rxn is the reverse of GSM rxn, 1 otherwise) */
$include "%gsm_rxn_pairs_path%"
v_exp_lb(gsm_j)
$include "%v_exp_lb_path%"
v_exp_ub(gsm_j)
$include "%v_exp_ub_path%"
kapp(j)
$include "../pPFBA_kapps_RBA.txt"
;

Variables
prosynSlackSum, inactiveFluxSum, fluxSum, z, v(j), venzSlack, fluxSlack, s_v_exp_lb(gsm_j), s_v_exp_ub(gsm_j), prosynSlackLB(pro), prosynSlackUB(pro)
;
venzSlack.lo = 0; venzSlack.up = %venzSlackAllow%;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
* 2e3 to allow changes in either direction
s_v_exp_lb.lo(gsm_j) = 0; s_v_exp_lb.up(gsm_j) = 2 * %vmax% * %nscale%;
s_v_exp_ub.lo(gsm_j) = 0; s_v_exp_ub.up(gsm_j) = 2 * %vmax% * %nscale%;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = %vmax% * %nscale%;
* bounds from GSM model
$include %gms_path%GSM_rxn_bounds.txt

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
v.fx(j)$prowaste(j) = 0; v.up('PROWASTE-TOTALPROTEIN') = inf; v.up('PROWASTE-PROTCYT') = inf;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = %vmax% * %nscale%;

$include "%phenotype_path%"

*** EQUATION DEFINITIONS ***
*Obj, Stoic, fluxSlackBounds, RiboCapacityNuc, RiboCapacityMito, UnknownRiboCapacity, Nonmodel
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, UnknownRiboCapacity
$include %enz_cap_declares_path%
;

* PSS..				prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
* Obj2..				inactiveFluxSum =e= sum(j$rxns_with_no_prodata(j), v(j));
Obj..				z =e= v('BIOSYN-PROTTOBIO');
* Obj3..				z =e= sum(j$rxns_metab(j), v(j));
*fluxSlackBounds..			fluxSlack =e= sum(gsm_j, s_v_exp_lb(gsm_j) + s_v_exp_ub(gsm_j));
Stoic(i)..			sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 		v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 		v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
UnknownRiboCapacity..	v('RIBOSYN-ribonuc') * %kribonuc% + v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * (sum(j$nuc_translation(j), NAA(j) * v(j)) + sum(j$mito_translation(j), NAA(j) * v(j)) + sum(j$unknown_ribo_translation(j), NAA(j) * v(j)));
*ProData(j)$prodata_set(j)..	v(j) =e= pro_val(j) * (1 - venzSlack);
* $include "../prosyn_abundance_constraints.txt"
* Inactive(j)$rxns_with_no_prodata(j)..	v(j) =e= 0;
NonModelProtAllo..	v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
* GSM upper and lower bounds for fluxes (if data available); slacks included in case necessary
* GSM_LB_exp(gsm_j)$v_exp_lb(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =g= (v_exp_lb(gsm_j) * %nscale%) - s_v_exp_lb(gsm_j);
* GSM_UB_exp(gsm_j)$v_exp_ub(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =l= (v_exp_ub(gsm_j) * %nscale%) + s_v_exp_ub(gsm_j);

* enforce measured kapps
*$include "../kapp_test.txt"
$include %enz_cap_eqns_path%
* $include "%gms_path%RBA_enzCapacityConstraints_declares_and_eqns_equality_version.txt"

*** BUILD OPTIMIZATION MODEL ***
Model rba
/all
/;
rba.optfile = 1;
* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Solve rba using lp minimizing z;

*** SOLVE ***
*Solve rba using lp minimizing fluxSlack;
*fluxSlack.up = fluxSlack.l + %epsilon%;
*Solve rba using lp minimizing z;

file ff /%system.FN%.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /%system.FN%.flux_gamsscaled.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt %vmin%),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
	);
);
putclose ff2;

file ff3 /%system.FN%.enzsyn.txt/;
ff3.nr = 2; put ff3;
loop(j,
	if ((rxns_enzsyn(j)),
		put j.tl:0, system.tab, (v.l(j)/%nscale%):0:15/;
	);
);

file ff3a /%system.FN%.flux_unscaled.txt/;
ff3a.nr = 2; put ff3a;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose ff3a;

