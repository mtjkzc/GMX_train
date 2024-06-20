#!/bin/bash

# Check if exactly one argument (file name) is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <cif_file>"
    exit 1
fi

# Get the filename from the argument
file="$1"

# Check if the file exists
if [ ! -f "$file" ]; then
    echo "Error: File $file does not exist."
    exit 1
fi

# Output file
output_file="${file%.cif}_CLEANED.cif"

# Temporary file to store filtered lines
temp_file=$(mktemp)

# Associative array to track the first occurrence of each condition
declare -A conditions

# Read the file line by line
while IFS= read -r line; do
    if [[ "$line" == *"ATOM "* && "$line" == *" A "* && "$line" == *"O OP3 "* ]]; then
        if [ -z "${conditions['A_OP3']}" ]; then
            conditions['A_OP3']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" A "* && "$line" == *"P P "* ]]; then
        if [ -z "${conditions['A_P']}" ]; then
            conditions['A_P']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" A "* && "$line" == *"O OP1 "* ]]; then
        if [ -z "${conditions['A_OP1']}" ]; then
            conditions['A_OP1']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" A "* && "$line" == *"O OP2 "* ]]; then
        if [ -z "${conditions['A_OP2']}" ]; then
            conditions['A_OP2']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" B "* && "$line" == *"O OP3 "* ]]; then
        if [ -z "${conditions['B_OP3']}" ]; then
            conditions['B_OP3']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" B "* && "$line" == *"P P "* ]]; then
        if [ -z "${conditions['B_P']}" ]; then
            conditions['B_P']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" B "* && "$line" == *"O OP1 "* ]]; then
        if [ -z "${conditions['B_OP1']}" ]; then
            conditions['B_OP1']=1
            continue
        fi
    elif [[ "$line" == *"ATOM "* && "$line" == *" B "* && "$line" == *"O OP2 "* ]]; then
        if [ -z "${conditions['B_OP2']}" ]; then
            conditions['B_OP2']=1
            continue
        fi
    fi
    # Write the line to the temporary file
    echo "$line" >> "$temp_file"
done < "$file"

# Move the temporary file to the output file
mv "$temp_file" "$output_file"

echo "Processed file saved as $output_file"

