#!/bin/bash

# Create the folder for storing ICTV Metadata
parent_dir=$(dirname "$PWD")
echo "$parent_dir/phallett"
mkdir -p "$parent_dir/phallett/data/Virus_Metadata_Resource"

#Delete existing "current" file 
if [ -f "$parent_dir/phallett/data/Virus_Metadata_Resource/current" ]; then
    rm "$parent_dir/phallett/data/Virus_Metadata_Resource/current"
fi

#Download latest current file 
cd "$parent_dir/phallett/data/Virus_Metadata_Resource"
curl -LO https://ictv.global/vmr/current

#Convert the file to a csv 
python3 <<EOF
import pandas as pd
df = pd.read_excel('current')
df.to_csv('VMR.csv', index=False)
EOF
