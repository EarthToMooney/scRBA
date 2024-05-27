************************* Run RBA model ********************
*       Author: Hoang Dinh
************************************************************

$INLINECOM /*  */
$include "./runRBA_GAMS_settings.txt"
* Scale values of all variables by a factor, then when write to file, descale them
$setGlobal nscale 1e3
$setGlobal nscaleback 1e-3

* For second run when the protein capacity is corrected, overwrite this default value of 1
* Value of 1 means the protein capacity is used freely, which leads to cheap proteins that
* catalyzed flux cycle to be produced.
* Flux cycle consists of the forward and reverse directions of a reversible reaction
* catalyzed by a protein that is cheaper than the dummy protein. Why is it cheaper?
* It is because the composition of that protein has higher fractions of cheap
* amino acids compared to the dummy protein which has the composition of the
* experimentally observations.
* This happens to system where protein availability
* is not the limiting factor (e.g., S. cerevisiae where nutrient or rRNA availability is
* the limiting factor).
$setGlobal corrected_protein_capacity_percentage 1

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 0.066

* Ribosome efficiency (3600 to convert /s into /h)
$setGlobal kribonuc 13.2*3600
$setGlobal kribomito 13.2*3600

* Enforce part of proteome allocate to non-modeled protein
$setGlobal nonmodeled_proteome_allocation 0.5 

options
	LP = cplex /*Solver selection*/
	limrow = 0 /*number of equations listed, 0 is suppresed*/
	limcol = 0 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 1000000 /*wall-clock time limit for solver in seconds*/
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
z, v(j), kapp_slack_ub(j),kapp_slack_lb(j)
;
kapp_slack_ub.lo(j) = 0;
kapp_slack_ub.up(j) = 1e5;
kapp_slack_lb.lo(j) = 0;
kapp_slack_lb.up(j) = 1e5;
kapp_slack_ub.fx(j)$(prosyn(j) or prowaste(j) or nuc_translation(j) or mito_translation(j) or uptake(j) or media(j)) = 0;
kapp_slack_lb.fx(j)$(prosyn(j) or prowaste(j) or nuc_translation(j) or mito_translation(j) or uptake(j) or media(j)) = 0;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e3 * %nscale%;

* Simulation top-level settings
* Enable or disable wasteful protein production, disabled by default (to solve faster)
* Note in solving: Enable protein waste flux might cause error for solver
* This is because enabling protein waste introduces several thousands more free variable to the system
* Thus, protein waste should only be implemented with actual data to constrain the free variable
v.fx(j)$prowaste(j) = 0;
*v.fx('PROSYN-rt7991') = 0.000001;
*v.fx('PROSYN-rt7221') = 0.000001;
*v.fx('PROSYN-rt3590') = 0.000001;
*v.fx('PROSYN-rt3599') = 0.000001;
*v.fx('PROSYN-rt5474') = 0.000001;
*v.fx('PROSYN-rt6544') = 0.000001;
*v.fx('PROSYN-rt5056') = 0.000001;
*v.fx('PROSYN-rtmATP6') = 0.000001;
*v.fx('PROSYN-rt2306') = 0.000001;
*v.fx('PROSYN-rt1572') = 0.000001;
*v.fx('PROSYN-rtmATP8') = 0.000001;
*v.fx('PROSYN-rtmATP9') = 0.000001;
*v.up('PROWASTE-rt7991') = 1000;
*v.lo('PROWASTE-rt7991') = 0;

* Disable all biomass reactions
* Condition-specific biomass reaction activation has to be done in phenotype.txt file
v.up('BIOSYN-BIODILAERO') = 0; 
v.up('BIOSYN-COFACTORANAEROBIC') = 0; v.up('BIOSYN-compCERANAEROBIC') = 0; v.up('BIOSYN-BIODILBATCHANAERO') = 0; v.up('BIOSYN-BIODILCHEMOANAERO') = 0; v.up('BIOSYN-BIODILSTARVE') = 0; v.up('BIOSYN-BIODILNOGAM') = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e3 * %nscale%;

$include "%phenotype_path%"
* Set all organism-specific parameters in phenotype.txt

* Set your NGAM in phenotype.txt since NGAM value depends on growth-condition
* Set your biomass composition in phenotype.txt since biomass composition depends on growth-condition

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, ModelProtAlloCorrection
;
*$include %enz_cap_declares_path%

Obj..			z =e= sum(j, kapp_slack_lb(j) + kapp_slack_ub(j));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
NonModelProtAllo..	v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
ModelProtAlloCorrection..	v('BIOSYN-PROTDUMMY2') =g= (1 - %corrected_protein_capacity_percentage%) * (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
*$include %kapp_slack_enz_cap_eqns_path%
$include "../../GAMS/model/kapp-slack-RBA_enzCapacityConstraints_eqns_equality_version.txt" 

*** BUILD OPTIMIZATION MODEL ***
Model rba
/all/;
*/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, NonModelProtAllo, MitoProtAllo, ModelProtAlloCorrection
*$include %enz_cap_declares_path%
*/;
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
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j) * %nscaleback%):0:15/;
	);
);
putclose ff2;

file ff3 /%system.FN%.kapp_info.txt/;
put ff3;
put 'rxn', system.tab, 'kapp_old', system.tab, 'kapp_slack_ub', system.tab, 'kapp_slack_lb'/;
loop(j,
	if ( (kapp_slack_ub.l(j) gt 0) or (kapp_slack_lb.l(j) gt 0),
		put j.tl:0, system.tab, kapp(j):0:15, system.tab, kapp_slack_ub.l(j):0:15, system.tab, kapp_slack_lb.l(j):0:15/;
	);
);
putclose ff3;

