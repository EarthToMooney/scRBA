* run RBA using default settings and outputs
$setGlobal gms %system.FN%
$include %model_root_path%GAMS/defaults.gms
* set up log files for printing results after runs
file bsr /'%model_root_path%binary_search_report.txt'/; put bsr;

for(max_iter=0 to %max_iter%,
	if(bi_UB - bi_LB > %bi_tol%,
		mu=(bi_UB+bi_LB)/2;
		solve rba using lp minimizing z;
		put 'mu =' mu 'status =';
		if(rba.modelstat=4,
			bi_UB=mu;
			put 'infeasible '/;
			putclose;
		else
			bi_LB=mu;
			put 'feasible (modelstat' rba.modelstat ')'/;
			putclose;
		);
	);
);
putclose;
$include %model_root_path%GAMS/default_outputs.gms
