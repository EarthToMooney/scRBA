* maximize ATPM to estimate GAM/NGAM
*** Minimize violation of fluxes that are zero due to non-production of enzymes ***
*       Author: Hoang Dinh
***********************************************************************************

$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 1000
$setGlobal venzSlackAllow 0
$setGlobal fluxSlackAllow 0
$setGlobal prosynSlackAllow 0

options
    LP = cplex /*Solver selection*/
    limrow = 1000000 /*number of equations listed, 0 is suppresed*/
    limcol = 1000000 /*number of variables listed, 0 is suppresed*/
    iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
    decimals = 8 /*decimal places for display statement*/
    reslim = 1000000 /*wall-clock time limit for solver in seconds*/
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
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
pro_val(j)
$include "%proteome_data_path%"
;

Variables
z, v(j), venzSlack(j), fluxslack(n), prosynSlackLB(pro), prosynSlackUB(pro)
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

* Disable enzyme synthesis and enzyme load network
v.fx(j)$rxns_enzsyn(j) = 0;
v.fx(j)$rxns_enzload(j) = 0;

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

* Carbohydrate, RNA, and protein fraction constraints
* C-lim protein fraction
*v.fx('BIOSYN-PROTTOBIO') = %mu% * (36.94 + 34.22*%mu%) / 100 * %nscale%;
* N-lim protein fraction
*v.fx('BIOSYN-PROTTOBIO') = %mu% * (10.57 + 108.56*%mu%) / 100 * %nscale%;
* RNA fraction for both C-lim and N-lim
*v.fx('BIOSYN-RNATOBIO') = %mu% * (3.79 + 13.88*%mu%) / 100 * %nscale%;

* Proteome allocation for purposes other than metabolism and ribosome
* Clim
*v.up('BIOSYN-PROTMODELED') = (1 - %nonmodeled_proteome_allocation%) * %mu% * (36.94 + 34.22*%mu%) / 100 * %nscale%;
* Nlim
*v.up('BIOSYN-PROTMODELED') = (1 - %nonmodeled_proteome_allocation%) * %mu% * (10.57 + 108.56*%mu%) / 100 * %nscale%;

* Additional constraints
*v.fx('RXN-EX_pyr_e_FWD-SPONT') = 0;
*v.fx('RXN-EX_btd_e_FWD-SPONT') = 0;

*** EQUATION DEFINITIONS ***
*$include "%fluxcap_declares_path%"
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, Nonmodel, MitoProtAllo
;

*Obj..              z =e= venzSlack;
*Obj..               z =e= sum(j$rxns_inactive(j), v(j)) + 1e6 * sum(j$prodata_set(j), venzSlack(j)) + sum(n,fluxSlack(n));
*Obj..               z =e= sum(j$rxns_inactive(j), v(j)) + 1e6 * sum(j$prodata_set(j), venzSlack(j));
v.up('RXN-ATPM_c_FWD-SPONT') = 1e4 * %nscale%;
v.lo('RXN-ATPM_c_FWD-SPONT') = 0;

Obj..               z =e= -v('RXN-ATPM_c_FWD-SPONT'); 
Stoic(i)..          sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito..      v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc..       v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
*ProData(j)$prodata_set(j).. v(j) =l= pro_val(j) * (1 + venzSlack(j)); 
*ProDataLB(j)$prodata_set(j)..   v(j) =g= pro_val(j) * (1 - venzSlack(j));
Nonmodel..          v('BIOSYN-PROTMODELED') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');

*$include "%fluxcap_path%"

*** BUILD OPTIMIZATION MODEL ***
Model rba
/all/;
*/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, ProData, ProDataLB, Nonmodel, MitoProtAllo
*$include "%fluxcap_declares_path%"
*/;
rba.optfile = 1;
*/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, ProData, Nonmodel

*** SOLVE ***
Solve rba using lp minimizing z

file ff /max_GAM.modelStat.txt/;
ff.nr = 2; put ff;
put rba.modelStat/;
putclose ff;

file ff2 /max_GAM.flux_gamsscaled.txt/;
ff2.nr = 2; put ff2;
loop(j,
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff2;

file ff3 /max_GAM.flux_essential_inactive_rxns_gamsscaled.txt/;
ff3.nr = 2; put ff3;
loop(j$rxns_inactive(j),
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff3;

file ff4 /max_GAM.venzSlack.txt/;
ff4.nr = 2; put ff4;
loop(j$prodata_set(j),
    if ( (venzSlack.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'venzSlack', system.tab, venzSlack.l(j):0:15/;
    );
);
putclose ff4;

file ff5 /max_GAM.fluxSlack.txt/;
ff5.nr = 2; put ff5;
loop(n,
    if ( (fluxSlack.l(n) gt 1e-12),
        put n.tl:0, system.tab, 'fluxSlack', system.tab, fluxSlack.l(n):0:15/;
    );
);
putclose ff5;

file ff6 /max_GAM.prosynSlack.txt/;
ff6.nr = 2; put ff6;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/; 
loop(pro,
    if ( prosynSlackUB.l(pro) > 1e-12 or prosynSlackLB.l(pro) > 1e-12,
        put pro.tl:0, system.tab, (100*prosynSlackUB.l(pro)-100*prosynSlackLB.l(pro)):0:15/;
    );
);
putclose ff6;
