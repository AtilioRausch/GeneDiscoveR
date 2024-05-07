#!/bin/bash

# Set the source and destination directories
input="path/to/results/OrthoFinder/*/Comparative_Genomics_Statistics/"
output="path/to/output/"

# Initial counter
counter=1

# Find all N0.tsv files in the source directory and its subdirectories
files=$(find $input -type f -name "Statistics_Overall.tsv")

# Iterate over the files and copy them to the destination directory with a unique name
for file in $files; do
  echo $file
  name=$(basename $file)
  newName="Statistics_Overall_$counter.tsv"
  cp $file $output$newName
  ((counter++))
done
