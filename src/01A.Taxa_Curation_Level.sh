#!/bin/bash

#Create a directory for taxa selected
parent_dir=$(dirname "$PWD")
mkdir -p "$parent_dir/phallett/data/Taxa_Selected"
#conda activate taxa_curation 
for genus_name in "$@";do
python3 "$parent_dir/phallett/src/01.Taxa_Curation_Level.py" $genus_name
done
