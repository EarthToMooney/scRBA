************************* Run RBA model ********************
*       Author: Hoang Dinh
************************************************************

$INLINECOM /*  */
$set path ./
$set mu 0.1

options
	solver = MINOS /*Solver selection*/
	limrow = 0 /*equation listing, 0 is suppresed*/
	limcol = 0 /*variable listing, 0 is suppresed*/
	optCR = 1e-15 /*relative optimality criterion solver default, important for problem with integer vars*/
	optCA = 1e-15 /*absolute optimality criterion solver default, important for problem with integer vars*/
	iterlim = 100000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 100000 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets

i
$include "%path%RBA_species.txt"
j
$include "%path%RBA_rxns.txt"
prosyn(j)
$include "%path%RBA_rxns_prosyn.txt"
nuctrnl(j)
$include "%path%RBA_rxns_nuctrnl.txt"
mitotrnl(j)
$include "%path%RBA_rxns_mitotrnl.txt"
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
v.lo(j) = 0;
v.up(j) = 1e4;

v.up('RXN-EX_glc__D_e_REV-SPONT') = 10;
v.up('RXN-EX_o2_e_REV-SPONT') = 1000;
v.fx('BIOSYN-BIODIL') = %mu%;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic
;

Obj..		z =e= sum(j$prosyn(j), v(j));
Stoic(i)..	sum(j, S(i,j)*v(j)) =e= 0;

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;
