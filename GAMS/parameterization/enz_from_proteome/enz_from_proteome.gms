******** Approximate enzyme concentration from proteome *********
*       Authors: (v1) Hoang Dinh, (v2) Eric Mooney
*****************************************************************

$INLINECOM /*  */
$include ../paths.txt
$include "./enz_from_proteome_GAMS_settings.txt"
$setGlobal nscale 1e6

* remove existing modelStat file to avoid issues w/ model status detection       
file ff /%system.FN%.modelStat.txt/; putclose ff '';

Sets
i
$include "%pro_and_enz_path%"
j
$include "%rxns_pro_and_enz_path%"
prodata(j)
$include "%rxns_pro_data_path%"
pronodata(j)
$include "%rxns_pro_nodata_path%"
enzout(j)
$include "%rxns_enz_path%"
dummy(j)
;

alias(j1,j)

Parameters
S(i,j)
$include "%sij_pro_and_enz_path%"
vprotexpmt(j)
$include "%proteome_data_path%"
;

Variables
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0;
v.up(j)$enzout(j) = 1e4;
v.up(j)$prodata(j) = vprotexpmt(j)*%nscale%;
v.up(j)$pronodata(j) = 0;

* protein abundance limits
$include "../protein_abundance_constraints.txt"

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic
;

Obj..			z =e= sum(j$dummy(j),v(j));
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;

*** BUILD OPTIMIZATION MODEL ***
Model optmodel
/all/;
optmodel.optfile = 1;

file ff3 /enz_flux_calculation.txt/;
ff3.nr=2; ff3.nz=1e-20;
put ff3;

* loop over each rxn, the objective to maximizing it and nothing else
loop(j1$(enzout(j1)),
    dummy(j)$(ord(j) ne ord(j1)) = no;
    dummy(j)$(ord(j) eq ord(j1)) = yes;
    display j1;
    
    Solve optmodel using lp maximizing z;
    If(optmodel.modelstat eq 1,
        put j1.tl:0 system.tab (v.l(j1)/%nscale%):0:15/;
	else 
		put j1.tl:0 system.tab 0/;
        );
	
    );

put ff;
put optmodel.modelStat/;
putclose;
