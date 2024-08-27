#!/bin/bash

# Replace '/path/to/your/directory' with the actual path
directory="."
master_directory="./pycore"
file_name="gsm_custom_functions.py"

# Find the first occurrence of the file as the master
master_file=$(find "$master_directory" -type f -name "$file_name" -print -quit)

if [ -z "$master_file" ]; then
  echo "No file found matching '$file_name'"
  exit 1
fi

echo "Master file: $master_file"

#find "$directory" -type f -name "$file_name" -exec bash -c '
#find "$directory" -type f -name "$file_name" ! -samefile "$master_file" -exec bash -c '
#  echo "Processing file: $0"
#  ln -f "$master_file" "$0"
#' {} \;
find "$directory" -type f -name "$file_name" ! -samefile "$master_file" -print0 | while IFS= read -r -d '' file; do
  echo "Processing file: $file"
  rm "$file"
  sleep 0.1
  ln -f "$master_file" "$file"
done
