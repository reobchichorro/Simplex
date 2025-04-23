#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p res

# Loop through all files in the 'files' directory
for file in files/*; do
    filename=$(basename "$file")           # Extract filename without path
    ./solve "$file" > "res/${filename}.txt"  # Run command and redirect output
done
