#!/bin/bash

# Get the name from the first argument
name="$1"

# Define the zip file and extraction directory
zip_file="fold_${name}.zip"
extraction_dir="${name}_af"

# Check if the zip file exists
if [ ! -f "$zip_file" ]; then
    echo "Error: File $zip_file does not exist."
    exit 1
fi

# Create the extraction directory if it doesn't exist
mkdir -p "$extraction_dir"

# Unzip the file into the extraction directory
unzip -q "$zip_file" -d "$extraction_dir"
mv $zip_file $extraction_dir

# Define the source file within the extracted directory
source_file="${extraction_dir}/fold_${name}_model_1.cif"
destination_file="./${name}_af_m0_PHOSPHATED.cif"
mv "$source_file" "$destination_file"
echo "= = = = = File moved and renamed to "$destination_file" = = = = ="

# Define the source file within the extracted directory
source_file="${extraction_dir}/fold_${name}_model_2.cif"
destination_file="./${name}_af_m1_PHOSPHATED.cif"
mv "$source_file" "$destination_file"
echo "= = = = = File moved and renamed to "$destination_file" = = = = ="

# Define the source file within the extracted directory
source_file="${extraction_dir}/fold_${name}_model_3.cif"
destination_file="./${name}_af_m2_PHOSPHATED.cif"
mv "$source_file" "$destination_file"
echo "= = = = = File moved and renamed to "$destination_file" = = = = ="

