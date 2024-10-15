*** Minimize violation of fluxes that are zero due to non-production of enzymes ***
*       Author: Hoang Dinh
***********************************************************************************

$INLINECOM /*  */
$include "./min_flux_sum_GAMS_settings.txt"
$setGlobal nscale 1000

options
	LP = soplex /*Solver selection*/
	limrow = 1000000 /*number of equations listed, 0 is suppresed*/
	limcol = 1000000 /*number of variables listed, 0 is suppresed*/
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
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxnsMetab(j) /*subset of reactions in the model that are metabolic*/
$include "%rxnsMetab_path%"
noEnzRxnSet(j) /*subset of index j of rxn that does not associate with enzyme*/
$include "%noEnzRxnSet_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
venzSyn(j)
$include "%enzSyn_venzSyn_path%"
;

Variables
z, v(j), venzSlack
;
venzSlack.lo = 0; venzSlack.up = %venzSlackAllow%;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4 * %nscale%;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e4 * %nscale%;

* Substrate and oxygenation, and secretions
$include "%phenotype_path%"

* Carbohydrate, RNA, and protein fraction constraints
v.fx('BIOSYN-CARBTOBIO') = %mu% * (0.4741 - 0.4662*%mu%) * %nscale%;
v.fx('BIOSYN-PROTTOBIO') = %mu% * (0.3736 + 0.3292*%mu%) * %nscale%;
*v.fx('BIOSYN-RNATOBIO') = %mu% * (0.0402 + 0.1370*%mu%) * %nscale%;

* Proteome allocation for purposes other than metabolism and ribosome
v.up('BIOSYN-PROTMODELED') = 0.55 * %mu% * (0.3736 + 0.3292*%mu%) * %nscale%;

* Limitation on proteome allocation to mitochondria
*v.up('BIOSYN-PROTMITO') = 0.5 * 0.55 * %mu% * (0.3736 + 0.3292*%mu%) * %nscale%;

* Additional constraints
*v.fx('RXN-EX_pyr_e_FWD-SPONT') = 0;
*v.fx('RXN-EX_btd_e_FWD-SPONT') = 0;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, Inactive
$include %enzSyn_declares_path%
;

Obj..			z =e= sum(j$rxnsMetab(j), v(j));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
Inactive..		sum(j$noEnzRxnSet(j), v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
$include %enzSyn_eqns_path%


*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, Inactive
$include %enzSyn_declares_path%
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

file ff /min_flux_sum.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /min_flux_sum.flux_gamsscaled.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 1e-12),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
	);
);
putclose ff2;
