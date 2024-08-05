#!/bin/bash
# for comparing yield prediction files to see if any settings have changed
output_max_withoutMP="./output_max_withoutMP"
output_max="./output_max"

comparison_method="sdiff"

exclude_files=("runRBA_max_prod.modelStat.txt" "report.txt" "RBA_result.json" "runRBA.flux.txt" "runRBA_max_prod.flux.txt" "flux.escher.csv")

function compare_directories() {
  dir1="$1"
  dir2="$2"

  for file in "$dir1"/*; do
    file2="${dir2}/${file##*/}"  # Get filename without path

    # Check if file is in the exclude list
    if [[ "$(echo "${exclude_files[@]}" | grep -w "${file##*/}")" ]]; then
      #echo "Skipping excluded file: $file"
      continue
    fi

    if [[ ! -f "$file2" ]]; then
      echo "File missing in $dir2: $file2"
      continue
    fi
    sdiff -sd "$file" "$file2" || echo "Files differ: $file and $file2"
    #"$comparison_method" "$file" "$file2" || echo "Files differ: $file and $file2"
    #diff -q "$file" "$file2"
  done
}

for dir in "$output_max_withoutMP"/*; do
  dir_name="${dir##*/}"  # Get directory name
  dir2="$output_max/$dir_name"

  if [[ ! -d "$dir2" ]]; then
    echo "Directory missing in output_max: $dir_name"
    continue
  fi

  compare_directories "$dir" "$dir2"
done

