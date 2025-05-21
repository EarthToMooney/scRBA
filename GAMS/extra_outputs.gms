* extra output files used after solving a model

file fi_v_sm / %gms%.SM_flux.txt /;
put fi_v_sm;
loop(sm_j,
	put sm_j.tl:0, system.tab, (v_sm.l(sm_j) / %nscale%):0:15/;
);