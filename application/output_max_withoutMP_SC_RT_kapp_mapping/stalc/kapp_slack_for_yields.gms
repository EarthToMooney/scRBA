* Run RBA model

$INLINECOM /*  */
$include "./runRBA_GAMS_settings.txt"
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1e5
$setGlobal vmax 1e3

* Optional constraint on allowed proteome allocation to mitochondrial proteins 
** (disable by setting to 1, enable in phenotype.txt file to ensure consistent behavior in subsequent models)
$setGlobal max_allowed_mito_proteome_allo_fraction 1

* Default value for ribosome efficiency
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

* if initial solve is infeasible, try adjustments (e.g., changing growth rate)
$setGlobal adjust_constraints_if_infeas 0 
* % of min experimental yield model must at least match
$setGlobal percent_min_exp_yield 0

* Enforce part of proteome allocate to non-modeled protein
** (disable by setting to 0, enable in phenotype.txt file to ensure consistent behavior in subsequent models)
$setGlobal nonmodeled_proteome_allocation 0
* default carbon slack (disabled by default)
$setGlobal carbonSlack 0

options
	LP = cplex /*Solver selection*/
	limrow = 0 /*number of equations listed, 0 is suppresed*/
	limcol = 0 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 600 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
prosyn(j)
$include "%prosyn_path%"
prowaste(j)
$include "%prowaste_path%"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_biomass(j)
$include "%biomass_path%"
rxns_add(j) /*list of heterologous rxns added for RBA model application*/
$include "%rxns_add_path%"
min_media(j)
$include "%model_root_path%GAMS/model/RBA_rxns_EXREV_minimal.txt"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
kapp(j)
$include "%kapp_path%"
mw(j)
$include "yield_denominators.txt"
mw_other_possible_substrates(j)
$include "yield_denominators_expanded.txt"
mu_current
;

Variables
z, prosynSlackSum, nonessentialInactiveFluxSum, inactiveFluxSum, fluxSum, v(j), venzSlack(j), fluxSlack, EnzLoadSlackPos(j), EnzLoadSlackNeg(j), slackSum
kapp_slack_ub(j),kapp_slack_lb(j)
mv_sub_sum "sum of (substrate flux * MW) for all substrates"
;
kapp_slack_ub.lo(j) = 0;
kapp_slack_ub.up(j) = inf;
kapp_slack_lb.lo(j) = 0;
kapp_slack_lb.up(j) = inf;
kapp_slack_ub.fx(j)$(prosyn(j) or prowaste(j) or nuc_translation(j) or mito_translation(j) or uptake(j) or media(j)) = 0;
kapp_slack_lb.fx(j)$(prosyn(j) or prowaste(j) or nuc_translation(j) or mito_translation(j) or uptake(j) or media(j)) = 0;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = %vmax% * %nscale%;
* bounds from GSM model
$include %model_root_path%GAMS/model/GSM_rxn_bounds.txt

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
v.fx(j)$prowaste(j) = 0; v.up('PROWASTE-TOTALPROTEIN') = inf; v.up('PROWASTE-PROTCYT') = inf;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;
* Turn off all heterologous rxns rxns_add
* Corresponding rxns for specific product will be turned on in phenotype.txt file
v.up(j)$rxns_add(j) = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = %vmax% * %nscale%;

$include "phenotype.txt"
* set mu as parameter for easier updating throughout the script
mu_current = %mu%;

* Set your NGAM in phenotype.txt since NGAM value depends on growth-condition
* Set your biomass composition in phenotype.txt since biomass composition depends on growth-condition

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, mvSubSum, minYield
$include %enz_cap_declares_path%
;

* Obj..			z =e= v('BIOSYN-PROTMODELED');
mvSubSum.. mv_sub_sum =e= sum(j$mw(j),v(j)*mw(j));
minYield.. minProductYield * mv_sub_sum * %percent_min_exp_yield% / 100 =l= %prod_mw% * v('%vprod%');
Obj..			z =e= sum(j, kapp_slack_lb(j) + kapp_slack_ub(j));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =e= mu_current * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =e= mu_current * sum(j$nuc_translation(j), NAA(j) * v(j));
NonModelProtAllo..	v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
$include "%gms_path%kapp-slack-RBA_enzCapacityConstraints_equality_version.txt" 

*** BUILD OPTIMIZATION MODEL ***
Model rba
/all
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

file ff /%system.FN%.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /%system.FN%.flux.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j) / %nscale%):0:15/;
	);
);
putclose ff2;

file ff3 /%system.FN%.kapp_info.txt/;
put ff3;
ff3.nr=2; ff3.nz=1e-30;
put 'rxn', system.tab, 'kapp_old', system.tab, 'kapp_slack_ub', system.tab, 'kapp_slack_lb'/;
loop(j,
	if ( (kapp_slack_ub.l(j) gt 0) or (kapp_slack_lb.l(j) gt 0),
		put j.tl:0, system.tab, kapp(j):0:15, system.tab, kapp_slack_ub.l(j):0:15, system.tab, kapp_slack_lb.l(j):0:15/;
	);
);
putclose ff3;

file ff4 /%system.FN%.kapp_slack_users.txt/;
put ff4;
loop(j,
	if ( (kapp_slack_ub.l(j) gt 0) or (kapp_slack_lb.l(j) gt 0),
		put "'", j.tl:0, "'"/;
	);
);
putclose ff4;
