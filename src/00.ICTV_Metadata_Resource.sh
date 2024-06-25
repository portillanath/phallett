#!/bin/bash

# Create the folder for storing ICTV Metadata
mkdir -p ~/phallett/data/Virus_Metadata_Resource

#Delete existing "current" file 
if [ -f "~/phallett/data/Virus_Metadata_Resource/current" ]; then
    rm ~/phallett/data/Virus_Metadata_Resource/current
fi

#Download latest current file 
cd ~/phallett/data/Virus_Metadata_Resource
wget --no-check-certificate https://ictv.global/vmr/current

#Convert the file to csv
python3 <<EOF
import pandas as pd
df = pd.read_excel('current')
df.to_csv('VMR.csv', index=False)
EOF
