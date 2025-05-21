* run RBA without bisection, using default settings and outputs
$if not setGlobal gms $setGlobal gms %system.FN%
$include %model_root_path%GAMS/defaults.gms
solve rba using lp minimizing z;
$include %model_root_path%GAMS/default_outputs.gms
