************************* Run RBA model ********************
*       Author: Hoang Dinh
************************************************************

$INLINECOM /*  */
$set path ../model/
$set mu 0.35
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
slack_apply(j)
$include "%path%RBA_rxns_slack_apply.txt"
;

Parameters
S(i,j)
$include "../model_zeroGAM/RBA_sij.txt"
NAA(j)
$include "%path%RBA_proteinLength.txt"
MW(j)
$include "%path%RBA_slack_species_MW.txt"
;

Variables
z, v(j), slL(j), slU(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4;
slL.fx(j) = 0;
slU.fx(j) = 0;
slL.lo(j)$slack_apply(j) = 0; slL.up(j)$slack_apply(j) = 1e4;
slU.lo(j)$slack_apply(j) = 0; slU.up(j)$slack_apply(j) = 1e4;

v.up('RXN-EX_glc__D_e_REV-SPONT') = 1000;
v.fx('RXN-EX_glc__D_e_FWD-SPONT') = 0;
v.up('RXN-EX_o2_e_REV-SPONT') = 1000;
v.fx('BIOSYN-BIODIL') = %mu%;
v.fx('RXN-EX_pyr_e_FWD-SPONT') = 0;
v.fx('RXN-EX_btd_e_FWD-SPONT') = 0;

slU.fx('ENZSYN-ATPS_m_FWD-ATPSCPLX') = 0;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic
$include %path%RBA_capacityConstraints_slacked_declares.txt
;

*Obj..		z =e= sum(j$prosyn(j), v(j));
Obj..		z =e= sum(j$slack_apply(j), slL(j) + slU(j));
Stoic(i)..	sum(j, S(i,j)*v(j)) =e= 0;
$include %path%RBA_capacityConstraints_slacked_eqns.txt

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic
$include %path%RBA_capacityConstraints_slacked_declares.txt
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

*** WRITE TO TEXT FILE ***
file ff /slack_active.txt/;
put ff;
loop(j$slack_apply(j),
	if ( (slL.l(j) gt 1E-10),
		put j.tl:0, system.tab, 'slL', system.tab, slL.l(j):0:11/;
	);
	if ( (slU.l(j) gt 1E-10),
		put j.tl:0, system.tab, 'slU', system.tab, slU.l(j):0:11/;
	);
);
putclose ff;

file ff2 /flux.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 1E-10),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:11/;
	);
);
putclose ff2;
