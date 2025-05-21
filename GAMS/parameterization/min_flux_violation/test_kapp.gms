$setGlobal gms %system.FN%

$INLINECOM /*  */
$include ../paths.txt
$offInline

$include %model_root_path%GAMS/RBA_single_iter.gms
file fi_vtAA /%gms%.v_trans_AA.txt/; put fi_vtAA;
put '/'/;
loop(AA,
	put "'" AA.tl:0 "' " (v_trans.l(AA) / v_trans_AA_sum.l):0:10/;
);
put '/'/;
putclose;
