* enzyme allocation
$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 1e3
* max fluxes allowed, and min fluxes deemed significant enough to report
$setGlobal vmax 1e4
$setGlobal vmin 1e-12
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal venzSlackAllow 0
$setGlobal fluxSlackAllow 0
$setGlobal prosynSlackAllow 0
* update to match solver tolerance
$setGlobal tol 2e-9

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
rxns_enzload_used(j)
$include "%rxns_enzload_used_path%"
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
dir(gsm_j,j) /* lists GSM rxn, RBA rxn, and direction (-1 if RBA rxn is the reverse of GSM rxn, 1 otherwise) */
$include "%gsm_rxn_pairs_path%"
v_exp_lb(gsm_j)
$include "%v_exp_lb_path%"
v_exp_ub(gsm_j)
$include "%v_exp_ub_path%"
;

* slacks for allowing fluxes to deviate from measured values when necessary
Variables
z, prosynSlackSum, nonessentialInactiveFluxSum, inactiveFluxSum, fluxSum, v(j), venzSlack(j), fluxSlack, s_v_exp_lb(gsm_j), s_v_exp_ub(gsm_j), prosynSlackLB(pro), prosynSlackUB(pro), EnzLoadSlackPos(j), EnzLoadSlackNeg(j), slackSum
;
venzSlack.lo(j) = 0; venzSlack.up(j) = %venzSlackAllow%;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
* 2e3 to allow changes in either direction
s_v_exp_lb.lo(gsm_j) = 0; s_v_exp_lb.up(gsm_j) = 2 * %vmax% * %nscale%;
s_v_exp_ub.lo(gsm_j) = 0; s_v_exp_ub.up(gsm_j) = 2 * %vmax% * %nscale%;
* slacks for enzyme load, to allow for more equal use of enzymes
EnzLoadSlackPos.lo(j) = 0; EnzLoadSlackPos.up(j) = %vmax% * %nscale%;
EnzLoadSlackNeg.lo(j) = 0; EnzLoadSlackNeg.up(j) = %vmax% * %nscale%;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = %vmax% * %nscale%;

** Disable enzyme load network for reactions that can't be used
v.fx(j)$rxns_enzload(j) = 0;
v.lo(j)$rxns_enzload_used(j) = %vmin%; 
v.up(j)$rxns_enzload_used(j) = %vmax% * %nscale%;

* Optional constraint on allowed proteome allocation to mitochondrial proteins (disable by setting to 1)
$setGlobal max_allowed_mito_proteome_allo_fraction 1
$setGlobal nonmodeled_proteome_allocation 0

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
Obj, PSS, Obj2, Obj3, Obj4, Stoic, RiboCapacityNuc, RiboCapacityMito, Nonmodel, GSM_LB, GSM_UB, fluxSlackBounds, MitoProtAllo
;
*$include "%fluxcap_declares_path%"

Obj..               z =e= v('PROWASTE-TOTALPROTEIN');
PSS..               prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
Obj2..				nonessentialInactiveFluxSum =e= sum(j$rxns_nonessential_inactive(j), v(j));
Obj3..				inactiveFluxSum =e= sum(j$rxns_inactive(j), v(j));
Obj4..				fluxSum =e= sum(j$rxns_metab(j), v(j));
Stoic(i)..          sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito..      v('RIBOSYN-ribomito') * %kribomito% =e= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc..       v('RIBOSYN-ribonuc') * %kribonuc% =e= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
Nonmodel..          v('BIOSYN-PROTMODELED') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..      v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
* GSM upper and lower bounds for fluxes (if data available); slacks included in case necessary
GSM_LB(gsm_j)$v_exp_lb(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =g= (v_exp_lb(gsm_j) * %nscale%) - s_v_exp_lb(gsm_j);
GSM_UB(gsm_j)$v_exp_ub(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =l= (v_exp_ub(gsm_j) * %nscale%) + s_v_exp_ub(gsm_j);
fluxSlackBounds..			fluxSlack =e= sum(gsm_j, s_v_exp_lb(gsm_j) + s_v_exp_ub(gsm_j));
*$include "%fluxcap_path%"

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Model prosyn /all/;
prosyn.optfile = 1;
Solve prosyn using lp minimizing prosynSlackSum;

prosynSlackSum.up = prosynSlackSum.l+%tol%;
* minimize disagreements with flux data
Model minFluxDeviations /all/;
minFluxDeviations.optfile = 1;
Solve minFluxDeviations using lp minimizing fluxSlack;

model ni_rxns /all/;
ni_rxns.optfile = 1;
Solve ni_rxns using lp minimizing nonessentialInactiveFluxSum;
v.up(j)$rxns_nonessential_inactive(j) = v.l(j)$rxns_nonessential_inactive(j)+%tol%;
*nonessentialInactiveFluxSum.up = nonessentialInactiveFluxSum.l+%tol%;
model i_rxns /all/;
i_rxns.optfile = 1;
Solve i_rxns using lp minimizing inactiveFluxSum;

*v.up(j)$rxns_inactive(j) = v.l(j)$rxns_inactive(j)+1e-8;
inactiveFluxSum.up = inactiveFluxSum.l+%tol%;
model m_rxns /all/;
m_rxns.optfile = 1;
Solve m_rxns using lp minimizing fluxSum;

* Solve again, encouraging more equal use of all enzymes
* force total flux to previous value
fluxSum.up = fluxSum.l+%tol%;
model min_prowaste /all/;
min_prowaste.optfile = 1;
Solve min_prowaste using lp minimizing z;


* apply enzyme load slacks
z.up = z.l+%tol%;
* Uncomment if using enzload values
$include "./enz_alloc_constraints.txt"
fluxSum.up = fluxSum.l+%tol%;
Equation Obj5; Obj5.. slackSum =e= sum(j, EnzLoadSlackNeg(j) + EnzLoadSlackPos(j));
Model eload /all/;
eload.optfile = 1;
Solve eload using lp minimizing slackSum;

file ff /enz_alloc.modelStat.txt/;
ff.nr = 2; put ff;
put eload.modelStat/;
putclose ff;

file ff2 /enz_alloc.flux_gamsscaled.txt/;
ff2.nr = 2; put ff2;
loop(j,
    if ( (v.l(j) gt %vmin%),
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

file ff6 /%system.FN%.s_v_exp.txt/;
ff6.nr = 2; ff6.pc=6; put ff6;
loop(gsm_j,
    if ( (s_v_exp_lb.l(gsm_j) gt %vmin%) or (s_v_exp_ub.l(gsm_j) gt %vmin%),
        put gsm_j.tl:0, s_v_exp_lb.l(gsm_j):0:15, s_v_exp_ub.l(gsm_j):0:15/;
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
