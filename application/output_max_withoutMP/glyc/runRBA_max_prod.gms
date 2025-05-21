* run RBA without bisection, using default settings and outputs
$setGlobal gms %system.FN%

$setGlobal settings_path "./runRBA_GAMS_settings.txt"

$include %model_root_path%GAMS/defaults.gms

* if initial solve is infeasible, try adjustments (e.g., changing growth rate)
$setGlobal adjust_constraints_if_infeas 0 
* % of min experimental yield model must at least match
$setGlobal percent_min_exp_yield 0
* default carbon slack (disabled by default)
$setGlobal carbonSlack 0

sets
rxns_add(j) /*list of heterologous rxns added for RBA model application*/
$include "%rxns_add_path%"
min_media(j)
$include "%gms_path%RBA_rxns_EXREV_minimal.txt"
;

parameters
mw(j)
$include "yield_denominators.txt"
mw_other_possible_substrates(j)
$include "yield_denominators_expanded.txt"
mu_current
;

variables
v_prod
mv_sub_sum "sum of (substrate flux * MW) for all substrates"
;

* Turn off all heterologous rxns rxns_add
* Corresponding rxns for specific product will be turned on in phenotype.txt file
v.up(j)$rxns_add(j) = 0;

equations
mvSubSum, minYield, maxProd
;
mvSubSum.. mv_sub_sum =e= sum(j$mw(j),v(j)*mw(j));
minYield.. minProductYield * mv_sub_sum * %percent_min_exp_yield% / 100 =l= %prod_mw% * v('%vprod%');
maxProd..			v_prod =e= -v('%vprod%');

Model rba
/all/;
rba.optfile = 1;

solve rba using lp minimizing v_prod;


* solve again, minimizing mv_sub_sum to find a better yield by discouraging wasteful substrate use 
* NOTE: (get rid of this if linear fractional programming version is developed)
* save current uptake bounds as parameters
Parameters
v_up(j), vprod_max
;
v_up(j) = v.up(j);
vprod_max = v.l('%vprod%');
v.lo('%vprod%') = vprod_max;
*if((rba.modelstat eq 1),
*	Solve rba using lp minimizing mv_sub_sum;
*);

* If not optimal, try again with bounds on vprod relaxed in case it's demanding too much.
* This isn't done initially in case setting specific bounds is important for replicating experimental results (e.g., for making multiple products).
v.up('%vprod%') = 1e3 * %nscale%;
v.lo('%vprod%') = 0 * %nscale%;

* If not optimal, find lowest uptakes that still allow optimal growth and exp. yields
* Retry while relaxing uptake upper bounds via slack variables, then minimize their sum
* Included in case uptake rate estimates are inaccurate.

* allow all possible uptakes in medium to reach max flux
v.up(j)$(uptake(j) and v_up(j) gt 0) = 1e3 * %nscale%;

Variables 
slack(j), slacksum
;
* default
slack.fx(j) = 0;
* allow slack for uptakes
slack.up(j)$(uptake(j) and v.up(j) gt 0) = inf;
slack.lo(j)$(uptake(j) and v.up(j) gt 0) = 0 * %nscale%;
* use initial uptake bounds as slack for other reactions
*$include "RBA_GAMS_defaults_from_FBA_initial.txt"
Equations slackcons, slackobj;
slackcons(j)$(slack.up(j) ne 0).. v(j) =l= v_up(j) + slack(j);
slackobj.. slacksum =e= sum(j, slack(j));
Model rba_slack /all/;

* added in case infeasibilities occur, and FBA predictions should be adjusted to match
Model fba /
Obj, Stoic, mvSubSum, minYield, slackcons, slackobj
/;

while((rba.modelstat ne 1 and (%adjust_constraints_if_infeas% ne 0)),
	rba_slack.optfile = 1;
	display 'UPDATE: find lowest uptakes that still allow optimal growth and exp. yields';
	Solve rba_slack using lp minimizing slacksum;
	if(rba_slack.modelstat eq 1, break;);

* If still not optimal, suggests minYield is wrong or the model's parameters don't perfectly match exp. conditions
* Try again, but find the highest yields under the original conditions
* turn off yield requirements
	minProductYield = 0;
	display 'UPDATE: find highest yields under the original conditions';
	Solve rba_slack using lp minimizing z;
	if(rba_slack.modelstat eq 1, break;);

* if still infeas., test if any media composition would work
	v.up(j)$(uptake(j)) = 1e3 * %nscale%;
	slack.up(j)$(uptake(j)) = 1e10 * %nscale%;
* find smallest uptake slacks needed for near-maximum flux
	display 'UPDATE: find smallest uptakes needed for near-maximum product flux';
	Solve rba_slack using lp minimizing slacksum;
	if(rba_slack.modelstat eq 1,
		v.fx("%vprod%") = v.l("%vprod%") * 0.99;
		Solve rba_slack using lp minimizing z;
	);
	break
);

$include %model_root_path%GAMS/default_outputs.gms


* write slacks and uptakes to file
file ff3 /runRBA.slack.txt/;
ff3.pc=6;
put ff3;
put 'rxn','flux','slack'/;
loop(j$(uptake(j) and v.up(j) gt 0),
	put j.tl:0, (v.l(j) / %nscale%):0:15,slack.l(j):0:15/;
);
putclose ff3;

if(rba.modelstat eq 1,
	abort.noError 'Optimal solution found';
);
if((%adjust_constraints_if_infeas% eq 0),
	abort.noError 'no optimal solution; set adjust_constraints_if_infeas to 1 to try tweaks'
);
* find smallest uptakes on minimal media
* set all other uptakes to 0
v.up(j)$uptake(j) = 0;
v.up(j)$(uptake(j) and min_media(j)) = 1e3 * %nscale%;
* set arbitrary uptake for glucose
*v.up('RXN-EX_glc__D_e_REV-SPONT') = 10 * %nscale%;

display 'UPDATE: find smallest uptakes needed for min exp. growth+yields on minimal media';
Solve rba_slack using lp minimizing slacksum;

v.lo('%vprod%') = 0;
v.up('%vprod%') = 1e3 * %nscale%;
* force slacks to not exceed RBA levels
slack.up(j) = slack.l(j);
* force uptakes to not exceed RBA levels
v.up(j)$uptake(j) = v.l(j);
display 'UPDATE: find FBA fluxes on minimal media using prior uptakes';
solve fba using lp minimizing z;
file ffba /runRBA.FBA_rerun_fluxes.txt/;
put ffba;
ffba.pc=6;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, 'v', (v.l(j) / %nscale%):0:15/;
	);
);
putclose ffba;

* output yield

set n "percent of growth rate from prior solutions" /10, 0/;
parameters
yield_denominator
;
file yields /runRBA.yields.txt/;
put yields;
yields.pc=6;
put 'prod', 'method', 'growth rate', 'yield'/;
yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
if(yield_denominator eq 0,
	yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
);
put '%vprod%','fba',(v.l('%biom_id%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
* find RBA yield under these conditions
solve rba using lp minimizing z;
yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
if(yield_denominator eq 0,
	yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
);
put '%vprod%','rba',(v.l('%biom_id%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
* allow PROWASTE rxns since without growth, there's no protein sink
v.up(j)$prowaste(j) = 1e3 * %nscale%;

loop(n,
	mu_current = mu * n.val/100;
	v.fx('%biom_id%') = mu_current * %nscale%;
	solve fba using lp minimizing z;
	yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
	if(yield_denominator eq 0,
		yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
	);
	put '%vprod%','fba',(v.l('%biom_id%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
	solve rba using lp minimizing z;
	yield_denominator = sum(j$mw(j),v.l(j)*mw(j));
* if original substrates are unused, assume glucose was used
	if(yield_denominator eq 0,
		yield_denominator = sum(j$mw_other_possible_substrates(j),v.l(j)*mw_other_possible_substrates(j));
	);
	put '%vprod%','rba',(v.l('%biom_id%') / %nscale%):0:15,(%prod_mw% * v.l('%vprod%') / yield_denominator):0:15/;
);

putclose yields;

file frba2 /runRBA.RBA_last_rerun_fluxes.txt/;
put frba2;
frba2.pc=6;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, 'v', (v.l(j) / %nscale%):0:15/;
	);
);
putclose frba2;
