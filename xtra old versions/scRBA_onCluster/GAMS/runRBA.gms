************************* Run RBA model ********************
*       Author: Hoang Dinh
************************************************************

$INLINECOM /*  */
$include "./runRBA_GAMS_settings.txt"

options
	LP = soplex /*Solver selection*/
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
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4;
v.fx('BIOSYN-BIODIL') = %mu%;

* Simulation top-level settings
* Enable or disable wasteful protein production, enabled by default
* v.fx(j)$prowaste(j) = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e4;

* Substrate and oxygenation
v.up('RXN-EX_glc__D_e_REV-SPONT') = 1e4;
v.fx('RXN-EX_glc__D_e_FWD-SPONT') = 0;

v.up('RXN-EX_o2_e_REV-SPONT') = 1e4;
v.up('RXN-EX_o2_e_FWD-SPONT') = 0;

* Comment out the line below if simulating anaerobic condition
* This is cofactor requirements excluding heme, which required oxygen to be synthesized
v.fx('BIOSYN-COFACTORANAEROBIC') = 0;
v.fx('BIOSYN-compCERANAEROBIC') = 0;
v.fx('BIOSYN-BIODILNOGAM') = 0;

* Carbohydrate, RNA, and protein fraction constraints
v.fx('BIOSYN-CARBTOBIO') = %mu% * (0.4741 - 0.4662*%mu%);
v.fx('BIOSYN-PROTTOBIO') = %mu% * (0.3736 + 0.3292*%mu%);
*v.fx('BIOSYN-RNATOBIO') = %mu% * (0.0402 + 0.1370*%mu%);

* Proteome allocation for purposes other than metabolism and ribosome
v.up('BIOSYN-PROTMODELED') = 0.55 * %mu% * (0.3736 + 0.3292*%mu%);

* Limitation on proteome allocation to mitochondria
*v.up('BIOSYN-PROTMITO') = 0.5 * 0.55 * %mu% * (0.3736 + 0.3292*%mu%);

* Additional constraints
*v.fx('RXN-EX_pyr_e_FWD-SPONT') = 0;
*v.fx('RXN-EX_btd_e_FWD-SPONT') = 0;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito
$include %enz_cap_declares_path%
;

Obj..			z =e= v('BIOSYN-PROTMODELED');
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 	v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 	v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
$include %enz_cap_eqns_path%

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito
$include %enz_cap_declares_path%
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

file ff /runRBA.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /runRBA.flux.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:11/;
	);
);
putclose ff2;
