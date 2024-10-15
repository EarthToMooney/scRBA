************************* Run RBA model ********************
*       Author: Hoang Dinh
************************************************************

$INLINECOM /*  */
$set path ./model/
$set mu 0.3
$set kribonuc 10.5*3600
$set kribomito 10.5*3600

options
	LP = soplex /*Solver selection*/
	limrow = 10000 /*number of equations listed, 0 is suppresed*/
	limcol = 0 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 1000000 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets

i
$include "%path%RBA_species.txt"
j
$include "%path%RBA_rxns.txt"
prosyn(j)
$include "%path%RBA_rxns_prosyn.txt"
nuc_translation(j)
$include "%path%RBA_nuc_translation.txt"
mito_translation(j)
$include "%path%RBA_mito_translation.txt"
uptake(j)
$include "%path%RBA_rxns_EXREV.txt"
media(j)
*$include "%path%RBA_rxns_EXREV_mineralMinimum.txt"
*$include "%path%RBA_rxns_EXREV_YNB.txt"
$include "%path%RBA_rxns_EXREV_YP.txt"
;

Parameters
S(i,j)
$include "%path%RBA_sij.txt"
NAA(j)
$include "%path%RBA_proteinLength.txt"
;

Variables
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4;
v.fx('BIOSYN-BIODIL') = %mu%;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e4;

* Substrate and oxygenation
v.up('RXN-EX_glc__D_e_REV-SPONT') = 1e4;
v.fx('RXN-EX_glc__D_e_FWD-SPONT') = 0;

v.up('RXN-EX_o2_e_REV-SPONT') = 1e4;
v.up('RXN-EX_o2_e_FWD-SPONT') = 0;

* Additional constraints
*v.fx('RXN-EX_pyr_e_FWD-SPONT') = 0;
*v.fx('RXN-EX_btd_e_FWD-SPONT') = 0;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic
$include %path%RBA_capacityConstraints_declares.txt
;

Obj..		z =e= sum(j$prosyn(j), v(j));
*Obj..		z =e= v('RXN-BIOMASS_SC_hvd_FWD-SPONT');
Stoic(i)..	sum(j, S(i,j)*v(j)) =e= 0;
$include %path%RBA_capacityConstraints_eqns.txt

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic
$include %path%RBA_capacityConstraints_declares.txt
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

file ff2 /flux.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 1E-10),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:11/;
	);
);
putclose ff2;
