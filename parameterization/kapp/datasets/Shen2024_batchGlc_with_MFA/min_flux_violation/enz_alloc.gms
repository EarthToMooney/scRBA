* finalize enzyme allocation

* define current gms file for use in other files
$setGlobal gms %system.FN%

* load default sets, parameters, and equations for kapp calculation
$include kapp_calc_defaults.gms

Sets
rxns_enzload_used(j)
$include "%rxns_enzload_used_path%"
rxns_NP_nonessential(j)
$include "./min_flux_violation.rxns_nonessential_with_no_prodata.txt"
;

Parameters
kapp(j)
$include "./kapp_init.txt"
prosynSlackSum_opt "optimal value of prosynSlackSum"
w_prowaste_opt
kappEstSlackSum_opt
fluxSlack_opt
fluxSum_j_NP_opt
fluxSum_opt
;

Variables
fluxSum_j_NP_nonEss
;

** Disable enzyme load network for reactions that can't be used
*v.fx(j)$rxns_enzload(j) = 0;
* force production of all enzymes involved in used rxns
v.lo(j)$rxns_enzload_used(j) = %vmin% * (1+1e-8); 
v.up(j)$rxns_enzload_used(j) = %vmax% * %nscale%;

Equations
flux_j_NP_NonEss, kappEst(j,j), kappMaxEst(j,j)
;

flux_j_NP_NonEss..				fluxSum_j_NP_nonEss =e= sum(j$rxns_NP_nonessential(j), v(j));
kappEst(j,j1)$enzload_rxn_coupling(j,j1).. (kappEstSlackPos(j1) - kappEstSlackNeg(j1) + v(j1)) * mu =e= v(j) * kapp(j1);
kappMaxEst(j,j1)$enzload_rxn_coupling(j,j1).. v(j1) * mu =l= v(j) * %kapp_max%;
$include "./enz_alloc_constraints.txt"

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Model kapp_calc /all-EnzCap-minPro-ProwasteLim/;
kapp_calc.optfile = 1;
Solve kapp_calc using lp minimizing prosynSlackSum;
put log; put 'minimized prosynSlackSum'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
prosynSlackSum_opt = prosynSlackSum.l;
prosynSlackSum.up = prosynSlackSum_opt;

Solve kapp_calc using lp minimizing w_prowaste;
put log; put 'minimized prowaste mass'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
w_prowaste_opt = w_prowaste.l;
w_prowaste.up = w_prowaste_opt;

Solve kapp_calc using lp minimizing kappEstSlackSum;
put log; put 'minimized kappEstSlackSum'/; putclose;
if (kapp_calc.modelStat ne 1, abort.noError "no optimal solution found";);
kappEstSlackSum_opt = kappEstSlackSum.l;
kappEstSlackSum.up = kappEstSlackSum_opt;

Model kapp_calc2 /all-EnzCap-minPro-ProwasteLim/;
kapp_calc2.optfile = 1;
Solve kapp_calc2 using lp minimizing fluxSlack;
put log; put 'minimized fluxSlack'/; putclose;
if (kapp_calc2.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSlack_opt = fluxSlack.l;
fluxSlack.up = fluxSlack_opt*(1+%tol%+.09);

Solve kapp_calc2 using lp minimizing fluxSum_j_NP;
put log; put 'minimized fluxSum_j_NP'/; putclose;
if (kapp_calc2.modelStat ne 1, abort.noError "no optimal solution found";);
fluxSum_j_NP_opt = fluxSum_j_NP.l;
* force rxns that were turned off to stay off
*v.fx(j)$(rxns_with_no_prodata(j) and (v.l(j) eq 0)) = 0;

fluxSum_j_NP.up = fluxSum_j_NP_opt*(1+%tol%);

*loop(enzload_rxn_coupling(j1,j),
*	if (v.l(j) eq 0, 
*			v.fx(j1) = 0;
*		);
*);

Solve kapp_calc2 using lp minimizing fluxSum;
put log; put 'minimized fluxSum'/; putclose;
if (kapp_calc2.modelStat ne 1, abort.noError "no optimal solution found";);
* force total flux to previous value
fluxSum.up = fluxSum.l*(1+%tol%);

* Solve again, encouraging more equal use of all enzymes
* NOTE: disabled this step since it can lead to arbitrary reductions in ENZLOAD fluxes, even when other ones aren't being used. This can increase kapps by reducing the denominator; how to fix this is unclear.
*Solve kapp_calc2 using lp minimizing slackSum;
*put log; put 'minimized uneven enzload distribution'/; putclose;

ff.nr = 2; put ff;
put kapp_calc2.modelStat/;
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

file enz_flux /%system.FN%.enz_flux_calculation.txt/;
enz_flux.nr = 2; put enz_flux;
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

