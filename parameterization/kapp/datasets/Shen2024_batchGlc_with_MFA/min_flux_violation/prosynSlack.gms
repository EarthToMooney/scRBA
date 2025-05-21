*** Minimize violation of fluxes that are zero due to non-production of enzymes ***
*       Authors: (v1) Hoang Dinh, (v2) Eric Mooney
***********************************************************************************

$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 1e5
* max fluxes allowed, and min fluxes deemed significant enough to report
$setGlobal vmax 1e3
$setGlobal vmin 1e-12
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal venzSlackAllow 0
$setGlobal prosynSlackAllow 0
* small value needed to ensure sequential problems aren't infeasible due to rounding errors
$setGlobal epsilon 1e-6

options
    LP = cplex /*Solver selection*/
    limrow = 1000000 /*number of equations listed, 0 is suppresed*/
    limcol = 1000000 /*number of variables listed, 0 is suppresed*/
    iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
    decimals = 8 /*decimal places for display statement*/
    reslim = 1000 /*wall-clock time limit for solver in seconds*/
    sysout = on /*solver status file report option*/
    solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
pro
$include "%unique_protein_set_path%"
gsm_j /* list of GSM model rxns */
$include "%gsm_rxns_path%"
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
unknown_ribo_translation(j)
$include "%unknown_ribo_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_biomass(j)
$include "%biomass_path%"
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
dir(gsm_j,j) /* lists GSM rxn, RBA rxn, and direction (-1 if RBA rxn is the reverse of GSM rxn, 1 otherwise) */
$include "%gsm_rxn_pairs_path%"
v_exp_lb(gsm_j)
$include "%v_exp_lb_path%"
v_exp_ub(gsm_j)
$include "%v_exp_ub_path%"
;

* slacks for allowing fluxes to deviate from measured values when necessary
Variables
prosynSlackSum, inactiveFluxSum, fluxSum, v(j), venzSlack(j), fluxSlack, s_v_exp_lb(gsm_j), s_v_exp_ub(gsm_j), prosynSlackLB(pro), prosynSlackUB(pro)
;
venzSlack.lo(j) = 0; venzSlack.up(j) = %venzSlackAllow%;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
* 2e3 to allow changes in either direction
s_v_exp_lb.lo(gsm_j) = 0; s_v_exp_lb.up(gsm_j) = 2 * %vmax% * %nscale%;
s_v_exp_ub.lo(gsm_j) = 0; s_v_exp_ub.up(gsm_j) = 2 * %vmax% * %nscale%;

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1
$setGlobal nonmodeled_proteome_allocation 0

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = %vmax% * %nscale%;
* bounds from GSM model
$include %gms_path%GSM_rxn_bounds.txt

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = %vmax% * %nscale%;

* Turning off all versions of biomass dilution reaction
* You need to turn on the respective version corresponding to your growth condition
v.fx(j)$rxns_biomass(j) = 0;
* protein abundance limits
$include "../prosyn_abundance_constraints.txt"

* Growth rate, substrate and oxygenation, and secretions
$include "%phenotype_path%"

*** EQUATION DEFINITIONS ***
Equations
Obj, Obj2, Obj3, Stoic, RiboCapacityNuc, RiboCapacityMito, UnknownRiboCapacity, Nonmodel, GSM_LB_exp, GSM_UB_exp, fluxSlackBounds, MitoProtAllo
;
Obj..				prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
Obj2..				inactiveFluxSum =e= sum(j$rxns_inactive(j), v(j));
Obj3..				fluxSum =e= sum(j$rxns_metab(j), v(j));
Stoic(i)..			sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityNuc..	v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
RiboCapacityMito..	v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
UnknownRiboCapacity..	v('RIBOSYN-ribonuc') * %kribonuc% + v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * (sum(j$nuc_translation(j), NAA(j) * v(j)) + sum(j$mito_translation(j), NAA(j) * v(j)) + sum(j$unknown_ribo_translation(j), NAA(j) * v(j)));
Nonmodel..			v('BIOSYN-PROTMODELED') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');

* GSM upper and lower bounds for fluxes (if data available); slacks included in case necessary
GSM_LB_exp(gsm_j)$v_exp_lb(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =g= (v_exp_lb(gsm_j) * %nscale%) - s_v_exp_lb(gsm_j);
GSM_UB_exp(gsm_j)$v_exp_ub(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =l= (v_exp_ub(gsm_j) * %nscale%) + s_v_exp_ub(gsm_j);
fluxSlackBounds..		fluxSlack =e= sum(gsm_j, s_v_exp_lb(gsm_j) + s_v_exp_ub(gsm_j));

* If any inactive essential rxns have kapps predicted from C1_calculate_kapp.py, 
* 	we will use them to constrain the fluxes so the effects of their protein demands are accounted for
*$include "../kapp_test.txt"
*$include "%fluxcap_path%"

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Model prosyn /all/;
prosyn.optfile = 1;
Solve prosyn using lp minimizing prosynSlackSum;

prosynSlackSum.up = prosynSlackSum.l + %epsilon%;

file ff /%system.FN%.modelStat.txt/;
ff.nr = 2; put ff; ff.pc=6;
put prosyn.modelStat/;
putclose ff;

file ff2 /%system.FN%.objectives.txt/;
ff2.nr = 2; put ff2; ff2.pc=6;
put 'prosynSlackSum',prosynSlackSum.l:0:11/;
put 'fluxSlackSum',(fluxSlack.l/%nscale%):0:11/;
put 'inactiveFluxSum',(inactiveFluxSum.l/%nscale%):0:11/;
put 'fluxSum',(fluxSum.l/%nscale%):0:11/;
putclose ff2;

file ff3 /%system.FN%.flux_gamsscaled.txt/;
ff3.nr = 2; put ff3;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff3;

file ff3a /%system.FN%.flux_unscaled.txt/;
ff3a.nr = 2; put ff3a;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose ff3a;

file ff4 /%system.FN%.flux_essential_inactive_rxns_gamsscaled.txt/;
ff4.nr = 2; put ff4;
loop(j$rxns_inactive(j),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose ff4;

file ff4a /%system.FN%.flux_essential_inactive_rxns_unscaled.txt/;
ff4a.nr = 2; put ff4a;
loop(j$rxns_inactive(j),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose ff4a;

file ff4b /%system.FN%.nonessential_inactive_rxns.txt/;
ff4b.nr = 2; put ff4b;
put '/'/;
loop(j$(rxns_inactive(j)),
    if ( (v.l(j) lt %vmin%),
        put j.tl:0/;
    );
);
put '/'/;
putclose ff4b;

file ff5 /%system.FN%.venzSlack.txt/;
ff5.nr = 2; put ff5;
loop(j$prodata_set(j),
    if ( (venzSlack.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'venzSlack', system.tab, venzSlack.l(j):0:15/;
    );
);
putclose ff5;

file ff6 /%system.FN%.s_v_exp.txt/;
ff6.nr = 2; ff6.pc=6; put ff6;
loop(gsm_j,
    if ( (s_v_exp_lb.l(gsm_j) gt 1e-12) or (s_v_exp_ub.l(gsm_j) gt 1e-12),
        put gsm_j.tl:0, s_v_exp_lb.l(gsm_j):0:15, s_v_exp_ub.l(gsm_j):0:15/;
    );
);
putclose ff6;

file ff7 /%system.FN%.prosynSlack.txt/;
ff7.nr = 2; put ff7;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/; 
loop(pro,
    if ( prosynSlackUB.l(pro) > %vmin% or prosynSlackLB.l(pro) > %vmin%,
        put pro.tl:0, system.tab, (100*prosynSlackUB.l(pro)-100*prosynSlackLB.l(pro)):0:15/;
    );
);
putclose ff7;
