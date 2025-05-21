* default output files used after solving a model; see defaults.gms for name details
put ff;
put rba.modelStat/;
putclose;

put fi_v;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j) / %nscale%):0:15/;
	);
);
putclose;

put fi_vgs;
loop(j,
	if ( (v.l(j) gt 0),
		put j.tl:0, system.tab, 'v', system.tab, (v.l(j)):0:15/;
	);
);
putclose;
