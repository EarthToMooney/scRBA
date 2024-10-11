* enzyme allocation

* define current gms file for use in other files
$setGlobal gms %system.FN%
* ignore constraints requiring production of measured but unused proteins from kapp calculations
$setGlobal ignore_measured_unused_constraints 1

$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 1e5
* max fluxes allowed, and min fluxes deemed significant enough to report
$setGlobal vmax 1e3
$setGlobal vmin 0 
* max predicted kapp, for use in approximating enzyme levels
$setGlobal kapp_max 1e30
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal fluxSlackAllow 0
$setGlobal prosynSlackAllow 0
* update to match solver tolerance
$setGlobal tol 3e-2

options
    LP = cplex /*Solver selection*/
    limrow = 1000000 /*number of equations listed, 0 is suppresed*/
    limcol = 1000000 /*number of variables listed, 0 is suppresed*/
    iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
    decimals = 8 /*decimal places for display statement*/
    reslim = 1000 /*wall-clock time limit for solver in seconds*/
    sysout = on /*solver status file report option*/
    solprint = on /*solution printing option*/
        
* remove existing modelStat file to avoid issues w/ model status detection       
file ff /%system.FN%.modelStat.txt/; putclose ff '';

Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
pro
$include "%unique_protein_set_path%"
prosyn(j)
$include "%prosyn_path%"
gsm_j /* list of GSM model rxns */
$include "%gsm_rxns_path%"
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
enzload_rxn_coupling(j,j)
$include "%gms_path%RBA_rxn_enzload_coupling.txt"
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
rxns_with_no_prodata(j) /*enzymatic rxns without proteomics data for all their enz subunits*/
$include "%rxns_with_no_prodata_path%"
enzload_with_no_prodata(j) /*enzload for enzymatic rxns without proteomics data for all their enz subunits*/
$include "./enzload_enz_with_no_prodata.txt"
rxns_NP_nonessential(j)
$include "./min_flux_violation.rxns_nonessential_with_no_prodata.txt"
prodata_set(j)
$include "%proteome_data_set_path%"
rxns_metab(j)
$include "%rxns_metab_path%"
;
alias (j1, j);

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
z, prosynSlackSum, kappEstSlackSum, fluxSum_j_NP_nonEss, fluxSum_j_NP, fluxSum, v(j), fluxSlack, s_v_exp_lb(gsm_j), s_v_exp_ub(gsm_j), prosynSlackLB(pro), prosynSlackUB(pro), EnzLoadSlackPos(j), EnzLoadSlackNeg(j), kappEstSlackPos(j), kappEstSlackNeg(j), slackSum
;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
* 2e3 to allow changes in either direction
s_v_exp_lb.lo(gsm_j) = 0; s_v_exp_lb.up(gsm_j) = 2 * %vmax% * %nscale%;
s_v_exp_ub.lo(gsm_j) = 0; s_v_exp_ub.up(gsm_j) = 2 * %vmax% * %nscale%;
* slacks for enzyme load, to allow for more equal use of enzymes
EnzLoadSlackPos.lo(j) = 0; EnzLoadSlackPos.up(j) = 10 * %vmax% * %nscale%;
EnzLoadSlackNeg.lo(j) = 0; EnzLoadSlackNeg.up(j) = 10 * %vmax% * %nscale%;
* slacks for apportioning enzyme load according to flux, to account for unmeasured enzymes that may not have been produced
kappEstSlackPos.lo(j) = 0; kappEstSlackPos.up(j) = inf;
kappEstSlackNeg.lo(j) = 0; kappEstSlackNeg.up(j) = inf;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = %vmax% * %nscale%;
* bounds from GSM model
$include %gms_path%GSM_rxn_bounds.txt

** Disable enzyme load network for reactions that can't be used
*v.fx(j)$rxns_enzload(j) = 0;
* force production of all enzymes involved in used rxns
v.lo(j)$rxns_enzload_used(j) = %vmin% * (1+1e-8); 
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
Obj, PSS, kappEstSlack, Obj2, Obj3, Obj4, Stoic, RiboCapacityNuc, RiboCapacityMito, UnknownRiboCapacity, Nonmodel, GSM_LB_exp, GSM_UB_exp, fluxSlackBounds, MitoProtAllo 
;

Obj..               z =e= v('PROWASTE-TOTALPROTEIN');
PSS..               prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
kappEstSlack..		kappEstSlackSum =e= sum(j, kappEstSlackPos(j) + kappEstSlackNeg(j));
Obj2..				fluxSum_j_NP_nonEss =e= sum(j$rxns_NP_nonessential(j), v(j));
Obj3..				fluxSum_j_NP =e= sum(j$rxns_with_no_prodata(j), v(j));
Obj4..				fluxSum =e= sum(j$rxns_metab(j), v(j));
Stoic(i)..          sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityNuc..       v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
RiboCapacityMito..      v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
UnknownRiboCapacity..	v('RIBOSYN-ribonuc') * %kribonuc% + v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * (sum(j$nuc_translation(j), NAA(j) * v(j)) + sum(j$mito_translation(j), NAA(j) * v(j)) + sum(j$unknown_ribo_translation(j), NAA(j) * v(j)));
Nonmodel..          v('BIOSYN-PROTMODELED') =e= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..      v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
* GSM upper and lower bounds for fluxes (if data available); slacks included in case necessary
GSM_LB_exp(gsm_j)$v_exp_lb(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =g= (v_exp_lb(gsm_j) * %nscale%) - s_v_exp_lb(gsm_j);
GSM_UB_exp(gsm_j)$v_exp_ub(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =l= (v_exp_ub(gsm_j) * %nscale%) + s_v_exp_ub(gsm_j);
fluxSlackBounds..		fluxSlack =e= sum(gsm_j, s_v_exp_lb(gsm_j) + s_v_exp_ub(gsm_j));
$include "./enz_alloc_constraints_with_kapp_estimates.txt"

$include "./enz_alloc_constraints.txt"
Equation Obj5; Obj5.. slackSum =e= sum(j, EnzLoadSlackNeg(j) + EnzLoadSlackPos(j));

file log /''/; 
* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Model rba /all/;
rba.optfile = 1;
Solve rba using lp minimizing prosynSlackSum;
put log; put 'minimized prosynSlackSum'/; putclose;
if (rba.modelStat ne 1, abort.noError "no optimal solution found";);
prosynSlackSum.up = prosynSlackSum.l*(1+(%tol%));

Solve rba using lp minimizing z;
put log; put 'minimized prowaste mass'/; putclose;
if (rba.modelStat ne 1, abort.noError "no optimal solution found";);
z.up = z.l*(1+(%tol%));

Solve rba using lp minimizing kappEstSlackSum;
put log; put 'minimized kappEstSlackSum'/; putclose;
if (rba.modelStat ne 1, abort.noError "no optimal solution found";);
kappEstSlackSum.up = kappEstSlackSum.l*(1+(1*%tol%));

Solve rba using lp minimizing fluxSlack;
put log; put 'minimized fluxSlack'/; putclose;
if (rba.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack.up = fluxSlack.l*(1+(1*%tol%));

Solve rba using lp minimizing fluxSum_j_NP;
put log; put 'minimized fluxSum_j_NP'/; putclose;
if (rba.modelStat ne 1, abort.noError "no optimal solution found";);
* force rxns that were turned off to stay off
v.fx(j)$(rxns_with_no_prodata(j) and (v.l(j) eq 0)) = 0;
fluxSum_j_NP.up = fluxSum_j_NP.l*(1+(1*%tol%));

*loop(enzload_rxn_coupling(j1,j),
*	if (v.l(j) eq 0, 
*			v.fx(j1) = 0;
*		);
*);

Solve rba using lp minimizing fluxSum;
put log; put 'minimized fluxSum'/; putclose;
if (rba.modelStat ne 1, abort.noError "no optimal solution found";);
* force total flux to previous value
fluxSum.up = fluxSum.l*(1+%tol%);

* Solve again, encouraging more equal use of all enzymes
* NOTE: disabled this step since it can lead to arbitrary reductions in ENZLOAD fluxes, even when other ones aren't being used. This can increase kapps by reducing the denominator; how to fix this is unclear.
*Solve rba using lp minimizing slackSum;
*put log; put 'minimized uneven enzload distribution'/; putclose;

ff.nr = 2; put ff;
put rba.modelStat/;
putclose;

file ff2 /%system.FN%.flux_gamsscaled.txt/;
ff2.nr = 2; put ff2;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose;

file ff2a /%system.FN%.flux_unscaled.txt/;
ff2a.nr = 2; put ff2a;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose;

file ff4a /%system.FN%.flux_essential_with_no_prodata_unscaled.txt/;
ff4a.nr = 2; put ff4a;
loop(j$rxns_with_no_prodata(j),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose;

file ff6 /%system.FN%.s_v_exp.txt/;
ff6.nr = 2; ff6.pc=6; put ff6;
loop(gsm_j,
    if ( (s_v_exp_lb.l(gsm_j) gt %vmin%) or (s_v_exp_ub.l(gsm_j) gt %vmin%),
        put gsm_j.tl:0, s_v_exp_lb.l(gsm_j):0:15, s_v_exp_ub.l(gsm_j):0:15/;
    );
);
putclose;

file ff7 /%system.FN%.enz_flux_calculation.txt/;
ff7.nr = 2; put ff7;
loop(j$(rxns_enzsyn(j) or rxns_enzload(j)),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, (v.l(j)/%nscale%):0:18/;
    );
);
putclose;

file ff8 /%system.FN%.enzload_used.txt/;
ff8.nr=2; put ff8;
loop(j$rxns_enzload_used(j),
    put j.tl:0, system.tab, v.l(j):0:15/;
);
putclose;

file ff8a /%system.FN%.kappEstSlack.txt/;
ff8a.nr=2; ff8a.pc=6; put ff8a;
put 'j','v','v_enzload','+ slack','- slack'/;
loop(enzload_rxn_coupling(j1,j),
	if ( (v.l(j) gt %vmin% and (kappEstSlackPos.l(j) gt 0 or kappEstSlackNeg.l(j) gt 0)),
		put j.tl:0, v.l(j):0:15, v.l(j1):0:15, kappEstSlackPos.l(j), kappEstSlackNeg.l(j)/;
	);
);
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff9 /%system.FN%.prosyn_unscaled.txt/;
ff9.nr=2; ff9.nz=1e-30; put ff9;
put '/'/;
loop(j$prosyn(j),
    put "'" j.tl:0 "' " (v.l(j)/%nscale%):0:15/;
);
put '/'/;
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff10 /%system.FN%.prosyn_gamsscaled.txt/;
ff10.nr=2; ff10.nz=1e-30; put ff10;
put '/'/;
loop(j$prosyn(j),
    put "'" j.tl:0 "' " v.l(j):0:15/;
);
put '/'/;
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff11 /%system.FN%.prosyn_frac_unscaled.txt/;
ff11.nr=2; ff11.nz=1e-30; put ff11;
put '/'/;
loop(j$prosyn(j),
    put "'" j.tl:0 "' " (v.l(j)/v.l('BIOSYN-PROTTOBIO')):0:15/;
);
put '/'/;
putclose;

* output protein levels for enforcing kapp calculation levels when needed
file ff12 /%system.FN%.prosyn_nonzero_gamsscaled.txt/;
ff12.nr=2; ff12.nz=1e-30; put ff12;
put '/'/;
loop(j$prosyn(j),
	if ( (v.l(j) ge %vmin%),	
		put "'" j.tl:0 "' " v.l(j):0:15/;
	);
);
put '/'/;
putclose;
