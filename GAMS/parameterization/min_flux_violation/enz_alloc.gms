* enzyme allocation
$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 100000
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal venzSlackAllow 0
$setGlobal fluxSlackAllow 0
$setGlobal prosynSlackAllow 0
* update to match solver tolerance
$setGlobal tol 1e-12

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
rxns_enzload_used(j)
$include "%rxns_enzload_used_path%"
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
rxns_nonessential_inactive(j)
$include "./min_flux_violation.nonessential_inactive_rxns.txt"
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
z, nonessentialInactiveFluxSum, inactiveFluxSum, fluxSum, v(j), venzSlack(j), fluxSlack(n), prosynSlackLB(pro), prosynSlackUB(pro), EnzLoadSlackPos(j), EnzLoadSlackNeg(j), slackSum
;
venzSlack.lo(j) = 0; venzSlack.up(j) = %venzSlackAllow%;
fluxSlack.lo(n) = 0; fluxSlack.up(n) = %fluxSlackAllow%;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
* slacks for enzyme load, to allow for more equal use of enzymes
EnzLoadSlackPos.lo(j) = 0; EnzLoadSlackPos.up(j) = 1e4 * %nscale%;
EnzLoadSlackNeg.lo(j) = 0; EnzLoadSlackNeg.up(j) = 1e4 * %nscale%;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4 * %nscale%;

** Disable enzyme load network for reactions that can't be used
v.fx(j)$rxns_enzload(j) = 0;
v.lo(j)$rxns_enzload_used(j) = %tol%; 
v.up(j)$rxns_enzload_used(j) = 1e4 * %nscale%;

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1
$setGlobal nonmodeled_proteome_allocation 0

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
Equations
Obj, Obj2, Obj3, Obj4, Stoic, RiboCapacityNuc, RiboCapacityMito, Nonmodel, MitoProtAllo
;
*$include "%fluxcap_declares_path%"

Obj..               z =e= v('PROWASTE-TOTALPROTEIN');
Obj2..				nonessentialInactiveFluxSum =e= sum(j$rxns_nonessential_inactive(j), v(j));
Obj3..				inactiveFluxSum =e= sum(j$rxns_inactive(j), v(j));
Obj4..				fluxSum =e= sum(j$rxns_metab(j), v(j));
Stoic(i)..          sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito..      v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc..       v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
Nonmodel..          v('BIOSYN-PROTMODELED') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..      v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');

*$include "%fluxcap_path%"

*** BUILD OPTIMIZATION MODEL ***
Model prosyn /all/;
prosyn.optfile = 1;

*** SOLVE ***
Solve prosyn using lp minimizing z;

* Solve again, encouraging more equal use of all enzymes
* force total flux to previous value
z.fx = z.l;

model ni_rxns /all/;
ni_rxns.optfile = 1;
Solve ni_rxns using lp minimizing nonessentialInactiveFluxSum;

nonessentialInactiveFluxSum.up = nonessentialInactiveFluxSum.l;
model i_rxns /all/;
i_rxns.optfile = 1;
Solve i_rxns using lp minimizing inactiveFluxSum;

inactiveFluxSum.up = inactiveFluxSum.l+%tol%;
model m_rxns /all/;
m_rxns.optfile = 1;
Solve m_rxns using lp minimizing fluxSum;

* Uncomment if using enzload values
$include "./enz_alloc_constraints.txt"
fluxSum.up = fluxSum.l+%tol%;
Equation Obj5; Obj5.. slackSum =e= sum(j, EnzLoadSlackNeg(j) + EnzLoadSlackPos(j));
Model eload /all/;
eload.optfile = 1;
Solve eload using lp minimizing slackSum;

file ff /enz_alloc.modelStat.txt/;
ff.nr = 2; put ff;
put m_rxns.modelStat/;
putclose ff;

file ff2 /enz_alloc.flux_gamsscaled.txt/;
ff2.nr = 2; put ff2;
loop(j,
    if ( (v.l(j) gt 1e-12),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff2;

file ff3 /enz_alloc.venzSlack.txt/;
ff3.nr = 2; put ff3;
loop(j$prodata_set(j),
    if ( (venzSlack.l(j) gt 0),
        put j.tl:0, system.tab, 'venzSlack', system.tab, venzSlack.l(j):0:15/;
    );
);
putclose ff3;

file ff5 /enz_alloc.fluxSlack.txt/;
ff5.nr = 2; put ff5;
loop(n,
    if ( (fluxSlack.l(n) gt 0),
        put n.tl:0, system.tab, 'fluxSlack', system.tab, fluxSlack.l(n):0:15/;
    );
);
putclose ff5;

file ff6 /enz_alloc.prosynSlack.txt/;
ff6.nr = 2; put ff6;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/;
loop(pro,
    if ( prosynSlackUB.l(pro) > 1e-12 or prosynSlackLB.l(pro) > 1e-12,
        put pro.tl:0, system.tab, (100*prosynSlackUB.l(pro)-100*prosynSlackLB.l(pro)):0:15/;
    );
);
putclose ff6;

file ff7 /enz_flux_calculation.txt/;
ff7.nr = 2;
ff7.nr = 2; put ff7;
loop(j$(rxns_enzsyn(j) or rxns_enzload(j)),
    if ( v.l(j) > 0,
        put j.tl:0, system.tab, v.l(j):0:15/;
    );
);
putclose ff7;

