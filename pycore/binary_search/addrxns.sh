#!/bin/bash

source_file="$1"
destination_file="$2"


contents=$(<"$source_file")
line_count=$(wc -l < "$destination_file")

insert_line=$((line_count - 1))

# Create a temporary file for filtered contents
temp_filtered=$(mktemp)

# Remove lines containing only '/' from the source file and save filtered contents
grep -v '^/$' "$source_file" > "$temp_filtered"

# Insert filtered contents into the second-to-last line of the destination file
sed -i "${insert_line}r $temp_filtered" "$destination_file"

# Remove duplicate lines from the destination file
sort -u -o "$destination_file" "$destination_file"

# Clean up temporary files
rm "$temp_filtered"

# Add a '/' back to the last line of the destination file
echo "/" >> "$destination_file"

echo "Contents of $source_file inserted and duplicates removed in the second-to-last line of $destination_file"
