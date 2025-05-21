* everything prior to model declaration for use in other files (e.g., min_flux_violation.gms, enz_alloc.gms)
* Created to help ensure consistency between kapp calculation steps
$INLINECOM /*  */

* ignore constraints requiring production of measured but unused proteins from kapp calculations
$setGlobal ignore_measured_unused_constraints 1
$setGlobal ignore_kapps 1
* max predicted kapp, for use in approximating enzyme levels
$setGlobal kapp_max 1e30

$include ../paths.txt
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal prosynSlackAllow 0
* small value needed to ensure sequential problems aren't infeasible due to rounding errors
$setGlobal epsilon 1e-6
* update to match solver tolerance
$setGlobal tol 1e-2

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
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
enzload_with_no_prodata(j) "enzload for enzymatic rxns without proteomics data for all their enz subunits"
$include "./enzload_enz_with_no_prodata.txt"
sm_j_lumped "rxns representing combinations of stoich. model rxns"
$include "../sm_j_lumped.txt"
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
sm_j_lumped_mappings(sm_j_lumped,sm_j) "grouping stoich. model rxns to their respective lumped counterparts"
$include "../sm_j_lumped_mappings.txt"
;

* slacks for allowing fluxes to deviate from measured values when necessary
Variables
w_prowaste, prosynSlackSum, kappEstSlackSum, fluxSum_j_NP, fluxSum, fluxSlack, s_v_exp_lb(sm_j), s_v_exp_ub(sm_j), s_v_exp_lb_lumped(sm_j_lumped), s_v_exp_ub_lumped(sm_j_lumped), prosynSlackLB(pro), prosynSlackUB(pro), EnzLoadSlackPos(j), EnzLoadSlackNeg(j), kappEstSlackPos(j), kappEstSlackNeg(j), slackSum
;
prosynSlackLB.lo(pro) = 0; prosynSlackLB.up(pro) = %prosynSlackAllow%;
prosynSlackUB.lo(pro) = 0; prosynSlackUB.up(pro) = %prosynSlackAllow%;
prosynSlackLB.up(pro)$unlimited_slack(pro) = 1;
prosynSlackUB.up(pro)$unlimited_slack(pro) = inf;
* 2e3 to allow changes in either direction
s_v_exp_lb.lo(sm_j) = 0; s_v_exp_lb.up(sm_j) = 2 * %vmax% * %nscale%;
s_v_exp_ub.lo(sm_j) = 0; s_v_exp_ub.up(sm_j) = 2 * %vmax% * %nscale%;
* slacks for lumped fluxes
s_v_exp_lb_lumped.lo(sm_j_lumped) = 0; s_v_exp_lb_lumped.up(sm_j_lumped) = 2 * %vmax% * %nscale%;
s_v_exp_ub_lumped.lo(sm_j_lumped) = 0; s_v_exp_ub_lumped.up(sm_j_lumped) = 2 * %vmax% * %nscale%;
* slacks for enzyme load, to allow for more equal use of enzymes
EnzLoadSlackPos.lo(j) = 0; EnzLoadSlackPos.up(j) = 10 * %vmax% * %nscale%;
EnzLoadSlackNeg.lo(j) = 0; EnzLoadSlackNeg.up(j) = 10 * %vmax% * %nscale%;
* slacks for apportioning enzyme load according to flux, to account for unmeasured enzymes that may not have been produced
kappEstSlackPos.lo(j) = 0; kappEstSlackPos.up(j) = inf;
kappEstSlackNeg.lo(j) = 0; kappEstSlackNeg.up(j) = inf;

Equations
Weight_Prowaste, kappEstSlack, ProsynSlacks, Obj1, Obj2, Obj3, sm_LB_exp, sm_UB_exp, fluxSlackBounds, lumped_flux(sm_j_lumped)
;

Weight_Prowaste..               w_prowaste =e= v('PROWASTE-TOTALPROTEIN');
kappEstSlack..		kappEstSlackSum =e= sum(j, kappEstSlackPos(j) + kappEstSlackNeg(j));
ProsynSlacks(pro)$(v_exp_pro(pro) gt 0).. sum(j$pro_prosyn(pro,j),v(j)) =e= v_exp_pro(pro) * (1 - prosynSlackLB(pro) + prosynSlackUB(pro)) * %nscale%;
Obj1..				prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
Obj2..				fluxSum_j_NP =e= sum(j$rxns_with_no_prodata(j), v(j));
Obj3..				fluxSum =e= sum(j$rxns_metab(j), v(j));
* sm upper and lower bounds for fluxes (if data available); slacks included in case necessary
sm_LB_exp(sm_j)$v_exp_lb(sm_j).. sum(j,dir(sm_j,j)*v(j)) =g= (v_exp_lb(sm_j) * %nscale%) - s_v_exp_lb(sm_j);
sm_UB_exp(sm_j)$v_exp_ub(sm_j).. sum(j,dir(sm_j,j)*v(j)) =l= (v_exp_ub(sm_j) * %nscale%) + s_v_exp_ub(sm_j);
fluxSlackBounds..		fluxSlack =e= sum(sm_j, s_v_exp_lb(sm_j) + s_v_exp_ub(sm_j)) + sum(sm_j_lumped, s_v_exp_lb_lumped(sm_j_lumped) + s_v_exp_ub_lumped(sm_j_lumped));
lumped_flux(sm_j_lumped).. sum(sm_j$sm_j_lumped_mappings(sm_j_lumped,sm_j),v_sm(sm_j)) =e= (smax(sm_j,sm_j_lumped_mappings(sm_j_lumped,sm_j)) * %nscale%) + s_v_exp_ub_lumped(sm_j_lumped) - s_v_exp_ub_lumped(sm_j_lumped);

Equation Obj5; Obj5.. slackSum =e= sum(j, EnzLoadSlackNeg(j) + EnzLoadSlackPos(j));

* allows command-line to output text; useful for debugging
file log /''/; 