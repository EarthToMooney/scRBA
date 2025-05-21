*** Minimize violation of fluxes that are zero due to non-production of enzymes ***
*       Authors: (v1) Hoang Dinh, (v2) Eric Mooney
***********************************************************************************
$INLINECOM /*  */

* define current gms file for use in other files
$setGlobal gms %system.FN%
* ignore constraints requiring production of measured but unused proteins from kapp calculations
$setGlobal ignore_measured_unused_constraints 1
$setGlobal ignore_kapps 1

$include ../paths.txt
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal prosynSlackAllow 0
* small value needed to ensure sequential problems aren't infeasible due to rounding errors
$setGlobal epsilon 1e-6 

$include "./min_flux_violation_GAMS_settings.txt"
* important to avoid errors when reading defaults file
$offInline

$include %model_root_path%GAMS/defaults.gms
v.up(j)$prowaste(j) = %vmax% * %nscale%;
Sets
rxns_with_no_prodata(j) "enzymatic rxns without proteomics data for all their enz subunits"
$include "%rxns_with_no_prodata_path%"
rxns_metab(j)
$include "%rxns_metab_path%"
unlimited_slack(pro) "proteins allowed to have unlimited slack."
;
*As in scRBA and rtRBA, allow slack to be used freely for ribosome subunits, to account for missing subunit paralogs or ambiguous stoichiometry.
unlimited_slack(pro)$sum(r,pro_ribo(r,pro)) = yes;

Parameters
v_exp_pro(pro)
$include "%v_exp_pro_path%"
v_exp_lb(sm_j)
$include "%v_exp_lb_path%"
v_exp_ub(sm_j)
$include "%v_exp_ub_path%"
;

* slacks for allowing fluxes to deviate from measured values when necessary
Variables
prosynSlackSum, fluxSum_j_NP, fluxSum, fluxSlack, s_v_exp_lb(sm_j), s_v_exp_ub(sm_j), prosynSlackLB(pro), prosynSlackUB(pro)
;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
prosynSlackLB.up(pro)$unlimited_slack(pro) = 0.3;
prosynSlackUB.up(pro)$unlimited_slack(pro) = inf;
* 2e3 to allow changes in either direction
s_v_exp_lb.lo(sm_j) = 0; s_v_exp_lb.up(sm_j) = 2 * %vmax% * %nscale%;
s_v_exp_ub.lo(sm_j) = 0; s_v_exp_ub.up(sm_j) = 2 * %vmax% * %nscale%;

Equations
ProsynSlacks, Obj1, Obj2, Obj3, sm_LB_exp, sm_UB_exp, fluxSlackBounds
;
ProsynSlacks(pro)$(v_exp_pro(pro) gt 0).. sum(j$pro_prosyn(pro,j),v(j)) =e= v_exp_pro(pro) * (1 - prosynSlackLB(pro) + prosynSlackUB(pro)) * %nscale%;
Obj1..				prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
Obj2..				fluxSum_j_NP =e= sum(j$rxns_with_no_prodata(j), v(j));
Obj3..				fluxSum =e= sum(j$rxns_metab(j), v(j));
* sm upper and lower bounds for fluxes (if data available); slacks included in case necessary
sm_LB_exp(sm_j)$v_exp_lb(sm_j).. sum(j,dir(sm_j,j)*v(j)) =g= (v_exp_lb(sm_j) * %nscale%) - s_v_exp_lb(sm_j);
sm_UB_exp(sm_j)$v_exp_ub(sm_j).. sum(j,dir(sm_j,j)*v(j)) =l= (v_exp_ub(sm_j) * %nscale%) + s_v_exp_ub(sm_j);
fluxSlackBounds..		fluxSlack =e= sum(sm_j, s_v_exp_lb(sm_j) + s_v_exp_ub(sm_j));

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
model minProSlack /all/;
minProSlack.optfile = 1;
Solve minProSlack using lp minimizing prosynSlackSum;
if (minProSlack.modelStat ne 1, abort.noError "no optimal solution found";);

prosynSlackSum.up = prosynSlackSum.l + %epsilon%;
* minimize disagreements with flux data
Model minFluxDeviations /all/;
minFluxDeviations.optfile = 1;
Solve minFluxDeviations using lp minimizing fluxSlack;
if (minFluxDeviations.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack.up = fluxSlack.l + (1e-7 + %epsilon%);

* minimize fluxes w/o proteomics data, to reduce reliance on rxns whose proteins aren't made
Model min_j_NP /all/;
min_j_NP.optfile = 1;
Solve min_j_NP using lp minimizing fluxSum_j_NP;
if (min_j_NP.modelStat ne 1, abort.noError "no optimal solution found";);
* force rxns that were turned off to stay off
v.fx(j)$(rxns_with_no_prodata(j) and (v.l(j) eq 0)) = 0;

fluxSum_j_NP.up = fluxSum_j_NP.l + (1e-4);

* minimize total flux sum, to satisfy parsimony assumption
Model minFlux /all/;
minFlux.optfile = 1;
Solve minFlux using lp minimizing fluxSum;

ff.nr = 2; put ff; ff.pc=6;
put minFlux.modelStat/;
putclose;

file ff2 /%system.FN%.objectives.txt/;
ff2.nr = 2; put ff2; ff2.pc=6;
put 'prosynSlackSum',prosynSlackSum.l:0:11/;
put 'fluxSlackSum',(fluxSlack.l/%nscale%):0:11/;
put 'fluxSum_j_NP',(fluxSum_j_NP.l/%nscale%):0:11/;
put 'fluxSum',(fluxSum.l/%nscale%):0:11/;
putclose;

file ff3 /%system.FN%.flux_gamsscaled.txt/;
ff3.nr = 2; put ff3;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
    );
);
putclose;

file ff3a /%system.FN%.flux_unscaled.txt/;
ff3a.nr = 2; put ff3a;
loop(j,
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose;

file ff4 /%system.FN%.flux_essential_with_no_prodata_gamsscaled.txt/;
file ff4c /%system.FN%.rxns_essential_with_no_prodata_gamsscaled.txt/;
ff4.nr = 2; put ff4;
ff4c.nr = 2; put ff4c;
put ff4c '/'/;
loop(j$rxns_with_no_prodata(j),
    if ( (v.l(j) gt %vmin%),
        put ff4 j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/ ff4c "'" j.tl:0 "'" system.tab v.l(j):0:15/;
    );
);
put ff4c '/';
putclose ff4 ff4c;

file ff4a /%system.FN%.flux_essential_with_no_prodata_unscaled.txt/;
ff4a.nr = 2; put ff4a;
loop(j$rxns_with_no_prodata(j),
    if ( (v.l(j) gt %vmin%),
        put j.tl:0, system.tab, 'v', system.tab, (v.l(j)/%nscale%):0:15/;
    );
);
putclose;

file ff4b /%system.FN%.rxns_nonessential_with_no_prodata.txt/;
ff4b.nr = 2; put ff4b;
put '/'/;
loop(j$(rxns_with_no_prodata(j)),
    if ( (v.l(j) le %vmin%),
        put j.tl:0/;
    );
);
put '/'/;
putclose;

file ff6 /%system.FN%.s_v_exp.txt/;
ff6.nr = 2; ff6.pc=6; put ff6;
loop(sm_j,
    if ( (s_v_exp_lb.l(sm_j) gt %vmin%) or (s_v_exp_ub.l(sm_j) gt %vmin%),
        put sm_j.tl:0, s_v_exp_lb.l(sm_j):0:15, s_v_exp_ub.l(sm_j):0:15/;
    );
);
putclose;

file ff7 /%system.FN%.prosynSlack.txt/;
ff7.nr = 2; put ff7;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/; 
loop(pro,
    if ( prosynSlackUB.l(pro) > %vmin% or prosynSlackLB.l(pro) > %vmin%,
        put pro.tl:0, system.tab, (100*prosynSlackUB.l(pro)-100*prosynSlackLB.l(pro)):0:15/;
    );
);
putclose;
