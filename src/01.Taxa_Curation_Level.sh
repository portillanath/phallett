#!/bin/bash

#Create a directory for taxa selected
mkdir -p ~/phallett/data/Taxa_Selected
#conda activate taxa_curation 
for genus_name in "$@";do
python3 ~/phallett/src/01.Taxa_Curation_Level.py $genus_name
done
#conda deactivate