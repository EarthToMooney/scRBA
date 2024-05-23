*** Minimize violation of fluxes that are zero due to non-production of enzymes ***
*       Author: Hoang Dinh
***********************************************************************************

$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 100000
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal venzSlackAllow 0
$setGlobal fluxSlackAllow 0
$setGlobal prosynSlackAllow 0

options
    LP = soplex /*Solver selection*/
    limrow = 1000000 /*number of equations listed, 0 is suppresed*/
    limcol = 1000000 /*number of variables listed, 0 is suppresed*/
    iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
    decimals = 8 /*decimal places for display statement*/
    reslim = 1000 /*wall-clock time limit for solver in seconds*/
    sysout = on /*solver status file report option*/
    solprint = on /*solution printing option*/
        
Sets
n /1*686/
pro
$include "%unique_protein_set_path%"
i
$include "%species_path%"
j
$include "%rxns_path%"
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_inactive(j)
$include "%rxns_inactive_path%"
prodata_set(j)
$include "%proteome_data_set_path%"
rxns_metab(j)
$include "%rxns_metab_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
pro_val(j)
$include "%proteome_data_path%"
;

* fluxSlack for allowing fluxes to deviate from measured values when necessary
Variables
prosynSlackSum, inactiveFluxSum, fluxSum, v(j), venzSlack(j), fluxSlack(n), prosynSlackLB(pro), prosynSlackUB(pro)
;
venzSlack.lo(j) = 0; venzSlack.up(j) = %venzSlackAllow%;
fluxSlack.lo(n) = 0; fluxSlack.up(n) = %fluxSlackAllow%;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1
$setGlobal nonmodeled_proteome_allocation 0

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4 * %nscale%;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e4 * %nscale%;

* Turning off all versions of biomass dilution reaction
* You need to turn on the respective version corresponding to your growth condition
v.fx('BIOSYN-COFACTORANAEROBIC') = 0; v.fx('BIOSYN-compCERANAEROBIC') = 0;
v.fx('BIOSYN-BIODILNOGAM') = 0; v.fx('BIOSYN-BIODILAERO') = 0; v.fx('BIOSYN-BIODILAERO-NOGAM') = 0; v.fx('BIOSYN-BIODILBATCHANAERO') = 0;
v.fx('BIOSYN-BIODILCHEMOANAERO') = 0; v.fx('BIOSYN-BIODILSTARVE') = 0;

* protein abundance limits
$include "../prosyn_abundance_constraints.txt"
* Growth rate, substrate and oxygenation, and secretions
$include "%phenotype_path%"

*** EQUATION DEFINITIONS ***
*$include "%fluxcap_declares_path%"
Equations
Obj, Obj2, Obj3, Stoic, RiboCapacityNuc, RiboCapacityMito, Nonmodel, MitoProtAllo
;
Obj..				prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
Obj2..				inactiveFluxSum =e= sum(j$rxns_inactive(j), v(j));
Obj3..				fluxSum =e= sum(j$rxns_metab(j), v(j));
Stoic(i)..			sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito..	v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc..	v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
Nonmodel..			v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');

* If any inactive essential rxns have kapps predicted from C1_calculate_kapp.py, 
* 	we will use them to constrain the fluxes so the effects of their protein demands are accounted for
*$include "../kapp_test.txt"
*$include "%fluxcap_path%"

* step 1: minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Model prosyn /all/;
prosyn.optfile = 1;
Solve prosyn using lp minimizing prosynSlackSum;

* step 2: minimize inactive fluxes, to reduce reliance on rxns whose proteins aren't made
prosynSlackSum.up = prosynSlackSum.l;
Model minInactive /all/;
minInactive.optfile = 1;
Solve minInactive using lp minimizing inactiveFluxSum;

* step 3: minimize total flux sum, to satisfy parsimony assumption
inactiveFluxSum.up = inactiveFluxSum.l;
Model minFlux /all/;
minFlux.optfile = 1;
Solve minFlux using lp minimizing fluxSum;

file ff /%system.FN%.modelStat.txt/;
ff.nr = 2; put ff;
put minFlux.modelStat/;
putclose ff;

file ff2 /%system.FN%.objectives.txt/;
ff2.nr = 2; put ff2;
put 'prosynSlackSum',system.tab,prosynSlackSum.l:0:11/;
put 'inactiveFluxSum',system.tab,(inactiveFluxSum.l/%nscale%):0:11/;
put 'fluxSum',system.tab,(fluxSum.l/%nscale%):0:11/;
putclose ff2;

file ff3 /%system.FN%.flux_gamsscaled.txt/;
ff3.nr = 2; put ff3;
loop(j,
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff3;

file ff4 /%system.FN%.flux_essential_inactive_rxns_gamsscaled.txt/;
ff4.nr = 2; put ff4;
loop(j$rxns_inactive(j),
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff4;

file ff3a /%system.FN%.flux_unscaled.txt/;
ff3a.nr = 2; put ff3a;
loop(j,
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose ff3a;

file ff4a /%system.FN%.flux_essential_inactive_rxns_unscaled.txt/;
ff4a.nr = 2; put ff4a;
loop(j$rxns_inactive(j),
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose ff4a;

file ff4b /%system.FN%.nonessential_inactive_rxns.txt/;
ff4b.nr = 2; put ff4b;
put '/'/;
loop(j$(rxns_inactive(j)),
    if ( (v.l(j) lt 1e-12),
        put j.tl:0/;
    );
);
put '/'/;
putclose ff4b;

file ff5 /%system.FN%.venzSlack.txt/;
ff5.nr = 2; put ff5;
loop(j$prodata_set(j),
    if ( (venzSlack.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'venzSlack', system.tab, venzSlack.l(j):0:15/;
    );
);
putclose ff5;

file ff6 /%system.FN%.fluxSlack.txt/;
ff6.nr = 2; put ff6;
loop(n,
    if ( (fluxSlack.l(n) gt 1e-12),
        put n.tl:0, system.tab, 'fluxSlack', system.tab, fluxSlack.l(n):0:15/;
    );
);
putclose ff6;

file ff7 /%system.FN%.prosynSlack.txt/;
ff7.nr = 2; put ff7;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/; 
loop(pro,
    if ( prosynSlackUB.l(pro) > 1e-12 or prosynSlackLB.l(pro) > 1e-12,
        put pro.tl:0, system.tab, (100*prosynSlackUB.l(pro)-100*prosynSlackLB.l(pro)):0:15/;
    );
);
putclose ff7;
