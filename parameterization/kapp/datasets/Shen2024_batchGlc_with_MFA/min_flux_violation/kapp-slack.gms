* Run RBA model

$INLINECOM /*  */
$include "../GAMS_setting_files/test_kapp_GAMS_settings.txt"
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1e5
$setGlobal vmax 1e3

* Optional constraint on allowed proteome allocation to mitochondrial proteins 
** (disable by setting to 1, enable in phenotype.txt file to ensure consistent behavior in subsequent models)
$setGlobal max_allowed_mito_proteome_allo_fraction 1

* Default value for ribosome efficiency
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

* Enforce part of proteome allocate to non-modeled protein
** (disable by setting to 0, enable in phenotype.txt file to ensure consistent behavior in subsequent models)
$setGlobal nonmodeled_proteome_allocation 0

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
pro
$include "%unique_protein_set_path%"
prosyn(j)
$include "%prosyn_path%"
prowaste(j)
$include "%prowaste_path%"
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
rxns_biomass(j)
$include "%biomass_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
kapp(j)
$include "%kapp_path%"
;

Variables
z, prosynSlackSum, nonessentialInactiveFluxSum, inactiveFluxSum, fluxSum, v(j), venzSlack(j), fluxSlack, prosynSlackLB(pro), prosynSlackUB(pro), EnzLoadSlackPos(j), EnzLoadSlackNeg(j), slackSum
kapp_slack_ub(j),kapp_slack_lb(j)
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
$include %gms_path%GSM_rxn_bounds.txt

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
*v.fx(j)$prowaste(j) = 0; v.up('PROWASTE-TOTALPROTEIN') = inf; v.up('PROWASTE-PROTCYT') = inf;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.fx(j)$rxns_biomass(j) = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = %vmax% * %nscale%;

$include "%phenotype_path%"
* Set all organism-specific and/or condition-specific parameters in phenotype.txt (e.g., NGAM and biomass composition)

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, UnknownRiboCapacity
$include %enz_cap_declares_path%
;

* Obj..			z =e= v('BIOSYN-PROTMODELED');
Obj..			z =e= sum(j, kapp_slack_lb(j) + kapp_slack_ub(j));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
UnknownRiboCapacity..	v('RIBOSYN-ribonuc') * %kribonuc% + v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * (sum(j$nuc_translation(j), NAA(j) * v(j)) + sum(j$mito_translation(j), NAA(j) * v(j)) + sum(j$unknown_ribo_translation(j), NAA(j) * v(j)));
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

