#!/bin/bash
# search for lines in RBA_sij.txt that contain a line from infes-prosyn.txt
infeas_rxns_path="infeas-rxns-3.txt"
#infeas_rxns_path="fva_all-infeas-names.txt"
feas_rxns_path="fva_all-feas.txt"
stoich_path="/storage/home/ejm6426/scratch/rtRBA-main/GAMS/model/RBA_sij.txt"
# grep -Ff $infeas_rxns_path $stoich_path
# remove all single quotes and everything after the first period
mets_used=$(grep -Ff $infeas_rxns_path $stoich_path | sed -E "s/'//g" | sed -E "s/\..*//" | sort -u)
# echo $mets_used

# # Determine contents of mets_used
# mets_used=$(awk 'NR==FNR { met=$1; sub(/\..*/, "", met); gsub("\047", "", met); mets[met]; next } { met=$1; sub(/\..*/, "", met); gsub("\047", "", met); if (!(met in mets)) print met }' infes-rxns.txt /Users/ejm6426/Desktop/rtRBA/rtRBA-main/GAMS/model/RBA_sij.txt | sort -u)

# # Print the resulting mets_used
# echo "$mets_used"

for met in $mets_used; do
    # Extracting text after the first period and before the first space, removing "'"
    # echo $met
    # search stoich_path for lines with $met, remove single quotes and everything before the first period
    matches=$(grep -w $met $stoich_path | sed -E 's/ .*//' |  sed -E "s/'//g" | sed -E "s/^[^\.]*\.//" | sort -u)
    # echo $matches

    # Checking if any of the matches is found in feas-rxns.txt
    found=false
    for match in $matches; do
        if grep -q "$match" "$feas_rxns_path"; then
            found=true
            break
        fi
    done

    # Printing the result based on the condition
    if [ "$found" == false ]; then
        echo "$met"
    fi
done

