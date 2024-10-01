#!/bin/bash
#for di in ls -d output_max/!(\(template\)); do rm -rf $di; done
dir=output_max_withoutMP
cd $dir
for di in $(ls -d *); do if [[ "$di" != "(template)" ]]; then rm -rf $di; fi; done
