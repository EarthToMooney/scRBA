* Find lowest ribosome elongation rate (by Eric Mooney)

* define current gms file for use in other files
$setGlobal gms %system.FN%
* ignore constraints requiring production of measured but unused proteins from kapp calculations
$setGlobal ignore_measured_unused_constraints 1

$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$include ../paths.txt
$include %model_root_path%GAMS/default_paths.txt
$setGlobal nscale 1e5
* max fluxes allowed, and min fluxes deemed significant enough to report
$setGlobal vmax 1e3
$setGlobal vmin 0
* slacks turned off by default, 
* 	but included to account for measurement errors when needed
$setGlobal prosynSlackAllow 0
* bisection tolerance
$setGlobal bi_tol .1

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
file log /''/;

Sets
t /* number of times bisection is tried */
/1*100/ 
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
r
$include %ribosomes_path%
j_ribo(r,j)
$include %j_ribo_path%
trans(r,j)
$include %translation_path%
trans_by_any_ribo(j)
$include %trans_by_any_ribo_path%
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
prodata_set(j)
$include "%proteome_data_set_path%"
rxns_metab(j)
$include "%rxns_metab_path%"
pro_ribo(r,pro)
$include %model_root_path%GAMS/model/pro_ribo_subunits.txt
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
kribo_bi
kribo(r)
$include %kribo_path%
iter_UB /* upper bound used for bisection */
iter_LB /* lower bound used for bisection */
mu
;
* highest known kribo_bi is 22 amino acids/ribosome/s, from https://doi.org/10.1128/ecosal.5.2.3
kribo_bi = 22; 
iter_LB = 0;
iter_UB = kribo_bi; /* default upper bound is highest kribo_bi */
iter_LB = 14.1;
iter_UB = 14.3; /* default upper bound is highest kribo_bi */

* slacks for allowing fluxes to deviate from measured values when necessary
Variables
prosynSlackSum, fluxSum_j_NP, fluxSum, v(j), fluxSlack, s_v_exp_lb(gsm_j), s_v_exp_ub(gsm_j), prosynSlackLB(pro), prosynSlackUB(pro),
v_ribo_total, v_prosyn_total
;
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
$include "../%phenotype_path%"
mu=%mu%;

*** EQUATION DEFINITIONS ***
Equations
Obj, Obj2, Obj3, Stoic, RiboCapacity, FlexibleRiboCapacity, Nonmodel, GSM_LB_exp, GSM_UB_exp, fluxSlackBounds, MitoProtAllo, prosynTotal, PSSlim
;
Obj..				prosynSlackSum =e= sum(pro, prosynSlackLB(pro) + prosynSlackUB(pro));
PSSlim..			sum(r,sum(pro$(not pro_ribo(r,pro)), prosynSlackLB(pro) + prosynSlackUB(pro))) =e= 0;
Obj2..				fluxSum_j_NP =e= sum(j$rxns_with_no_prodata(j), v(j));
Obj3..				fluxSum =e= sum(j$rxns_metab(j), v(j));
Stoic(i)..			sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacity(r).. sum(j$j_ribo(r,j),v(j) * kribo(r)) =g= mu * sum(j$trans(r,j), NAA(j) * v(j));
FlexibleRiboCapacity.. sum(r,sum(j$j_ribo(r,j),v(j) * kribo(r))) =g= mu * sum(r,sum(j$trans(r,j), NAA(j) * v(j))) + sum(j$trans_by_any_ribo(j), NAA(j) * v(j));
Nonmodel..			v('BIOSYN-PROTMODELED') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');
MitoProtAllo..		v('BIOSYN-PROTMITO') =l= %max_allowed_mito_proteome_allo_fraction% * v('BIOSYN-PROTMODELED');
prosynTotal..		v_prosyn_total =e= v('BIOSYN-PROTTOBIO');

* GSM upper and lower bounds for fluxes (if data available); slacks included in case necessary
GSM_LB_exp(gsm_j)$v_exp_lb(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =g= (v_exp_lb(gsm_j) * %nscale%) - s_v_exp_lb(gsm_j);
GSM_UB_exp(gsm_j)$v_exp_ub(gsm_j).. sum(j,dir(gsm_j,j)*v(j)) =l= (v_exp_ub(gsm_j) * %nscale%) + s_v_exp_ub(gsm_j);
fluxSlackBounds..		fluxSlack =e= sum(gsm_j, s_v_exp_lb(gsm_j) + s_v_exp_ub(gsm_j));

* minimize disagreement with proteomics data, while allowing some where needed (e.g., measurement errors)
Model minProSlack /all/;
minProSlack.optfile = 1;
Solve minProSlack using lp minimizing prosynSlackSum;
if (minProSlack.modelStat ne 1, abort.noError "no optimal solution found";);
loop(pro$(not pro_ribo(r,pro)),
	prosynSlackUB.up(pro) = prosynSlackUB.l(pro);
	prosynSlackLB.up(pro) = prosynSlackLB.l(pro);
);
* show the initial bounds
put log; put 'UB:' iter_UB ' LB:' iter_LB/; putclose;

*while the RBA model status is infeasible, or the model bounds are not sufficiently close
loop(t$(iter_UB - iter_LB > %bi_tol% or minProSlack.modelstat = 4),
  /*pick a point equidistant between the two bounds*/
	kribo_bi = (iter_UB + iter_LB)/2;
	kribo(r) = kribo_bi*3600;

	SOLVE minProSlack USING lp MINIMIZING prosynSlackSum;

  /*based on solution status update the bounds*/
	if ( minProSlack.modelstat = 1 ,
		iter_UB = kribo_bi;
	else
		iter_LB = kribo_bi;
	);
	put log; put 'UB:' iter_UB ' LB:' iter_LB/; putclose;
	display iter_UB, iter_LB;
)

put ff;
ff.nr = 2; put ff; ff.pc=6;
put minProSlack.modelStat/;
putclose;

file kr /"../kribo.txt"/; put kr;
put '/'/;
loop(r,
	put "'"r.tl:0"'" kribo(r)/;
);
put '/';
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

file ff7 /%system.FN%.prosynSlack.txt/;
ff7.nr = 2; put ff7;
put 'index', system.tab, 'prosynSlack(% higher/lower than measured value)'/; 
loop(pro,
    if ( prosynSlackUB.l(pro) > %vmin% or prosynSlackLB.l(pro) > %vmin%,
        put pro.tl:0, system.tab, (100*prosynSlackUB.l(pro)-100*prosynSlackLB.l(pro)):0:15/;
    );
);
putclose;