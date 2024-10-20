#!/bin/bash
#for di in ls -d output_max/!(\(template\)); do rm -rf $di; done
#dir=output_max_withoutMP
#cd $dir
file='runRBA.py'
# file='runRBA_GAMS_settings.txt'
for di in $(ls -d output_max*/*/$file); do
# for di in $(ls -d output_max_withoutMP_SC_RT_kapp_mapping*/*/$file); do
	if [[ "$di" != "(template)" && "$di" != "run.sh" ]]; then
		cp output_max_withoutMP/\(template\)/$file $di
		# cp output_max_withoutMP_SC_RT_kapp_mapping/\(template\)/$file $di
	fi
done