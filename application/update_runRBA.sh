#!/bin/bash
#for di in ls -d output_max/!(\(template\)); do rm -rf $di; done
#dir=output_max_withoutMP
#cd $dir
for di in $(ls -d output_max*/*); do
	if [[ "$di" != "(template)" && "$di" != "run.sh" ]]; then
		cp output_max_withoutMP/\(template\)/runRBA.py $di
	fi
done
