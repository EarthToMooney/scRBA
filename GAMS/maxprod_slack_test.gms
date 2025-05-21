************************* Run RBA model ********************
*       Authors: (v1) Hoang Dinh, (v2) Eric Mooney
************************************************************

$INLINECOM /*  */
$include "./runRBA_GAMS_settings.txt"
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1e5

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1

* Ribosome efficiency
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

* Enforce part of proteome allocate to non-modeled protein
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
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e3 * %nscale%;
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
v.up(j)$media(j) = 1e3 * %nscale%;

$include "phenotype.txt"

* Set your NGAM in phenotype.txt since NGAM value depends on growth-condition
* Set your biomass composition in phenotype.txt since biomass composition depends on growth-condition

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo
$include %enz_cap_declares_path%
;

Obj..			z =e= -v('%vprod%');
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
NonModelProtAllo..	v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
$include %enz_cap_eqns_path%

*** BUILD OPTIMIZATION MODEL ***
Model rba
/all/;
rba.optfile = 1;

*** SOLVE ***
*Solve rba using lp minimizing z;

* If not optimal, find lowest uptakes that still allow optimal growth and exp. yields
* Retry while relaxing uptake upper bounds via slack variables, then minimize their sum
* Included in case uptake rate estimates are inaccurate.

* save current uptake bounds as parameters
Parameter 
v_up(j)
;
v_up(j) = v.up(j);
* allow all possible uptakes in medium to reach max flux
v.up(j)$(uptake(j)) = 1e10 * %nscale%;

Variables 
slack(j), slacksum
;
* default
slack.fx(j) = 0;
* allow slack for uptakes
slack.up(j)$(uptake(j)) = 1e10 * %nscale%;
*slack.lo(j)$(uptake(j) and v.up(j) gt 0) = 0 * %nscale%;
* use initial uptake bounds as slack for other reactions
*$include "RBA_GAMS_defaults_from_FBA_initial.txt"
Equations slackcons, slackobj;
slackcons(j)$(slack.up(j) ne 0).. v(j) =l= v_up(j) + slack(j);
slackobj.. slacksum =e= sum(j, slack(j));
Model rba_slack /all/;

rba_slack.optfile = 1;
Solve rba_slack using lp minimizing slacksum;

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

* write slacks and uptakes to file
file ff3 /%system.FN%.slack.txt/;
put ff3;
loop(j$(uptake(j) and v.up(j) gt 0),
	put j.tl:0, system.tab, (v.l(j) / %nscale%):0:15, system.tab, 'slack', system.tab, slack.l(j):0:15/;
);
