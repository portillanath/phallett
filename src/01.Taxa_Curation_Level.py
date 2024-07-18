#!/usr/bin/env python
# coding: utf-8

# In[5]:

import os
import time
import pandas as pd
from Bio import Entrez
from pathlib import Path

def retrieve_genomes(genera_list):
    # Read the genus or bunch of genera to run
    current_dir = Path.cwd()
    parent_dir = current_dir.parent
    print(parent_dir)
    ICTV_assignation = pd.read_csv(f"{parent_dir}/phallett/data/Virus_Metadata_Resource/VMR.csv")
    ICTV_assignation = ICTV_assignation.rename(columns=lambda x: x.strip())  # Remove leading/trailing whitespaces
    ICTV_assignation.columns = ["Sort","Isolate_Sort","Realm","Subrealm","Kingdom","Subkingdom","Phylum","Subphylum","Class","Subclass","Order","Suborder","Family","Subfamily","Genus",
"Subgenus","Species","Exemplar_or_additional_isolate","Virus_name(s)","Virus_name_abbreviation(s)","Virus_isolate_designation","Virus_GENBANK_accession","Virus_REFSEQ_accession",
"Genome_coverage","Genome_composition","Host_source"]
    ICTV_assignation = ICTV_assignation[(ICTV_assignation["Host_source"] == "archaea") | (ICTV_assignation["Host_source"] == "bacteria")]
    ICTV_assignation = ICTV_assignation[(ICTV_assignation["Genome_coverage"] == "Complete genome") | (ICTV_assignation["Genome_coverage"] == "Complete coding genome")]
    ICTV_assignation.to_csv(f"{parent_dir}/phallett/test/ICTV_assignation_Complete.csv", index=False)

    # Extract accessions for search on NCBI
    accessions = ICTV_assignation["Virus_GENBANK_accession"].unique()
    print(f"The number of phages with complete genome for all ICTV database is: {len(accessions)}")

    # Creates a folder for storage the selected taxas
    import os

# Creates a folder for storage the selected taxas
    ncbi_genome_actual = f"{parent_dir}/phallett/data/Taxa_Selected"
    ncbi_genome_actual = os.path.expanduser(ncbi_genome_actual)
    os.makedirs(ncbi_genome_actual, exist_ok=True)

    for genus in genera_list:
        # Filter the ICTV assignation for the current genus
        genus_data = ICTV_assignation[ICTV_assignation["Genus"] == genus]

        # Check if there are any records found for the genus
        if genus_data.empty:
            print(f"No records found for the genus: {genus}")
            continue

        # Create a subfolder for the genus if it doesn't exist
        genus_folder = os.path.join(ncbi_genome_actual, genus)
        print(f"Trying to create folder:{genus_folder}")
        os.makedirs(genus_folder, exist_ok=True)

        # Iterate over accessions for the current genus
        for accession in genus_data["Virus_GENBANK_accession"]:
            file_name = f"{accession}.fasta"
            file_path = os.path.join(genus_folder, file_name)

            while True:
                try:
                    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
                        seq = handle.read()
                        if seq:
                            with open(file_path, 'w') as file:
                                file.write(seq)
                            print(f"Downloaded {file_name}")
                            break
                except Exception as e:
                    print(f"Error: {str(e)}")
                    time.sleep(5)
                    continue

                # Check if the file was created
                if os.path.exists(file_path):
                    print(f"New FASTA file created: {file_path}")
    
if __name__ == "__main__":
    import sys

    # Check if the command line arguments are provided
    if len(sys.argv) == 1:
        print("Please provide either a CSV file path with the list of genus or type space-separated genus names")
        sys.exit(1)

    # Get the command-line arguments
    input_argument = sys.argv[1]

    if os.path.exists(input_argument):
        # Read the text file with genus names
        with open(input_argument, 'r') as file:
            genera_list_txt = [line.strip() for line in file.readlines() if line.strip()]

        if not genera_list_txt:
            print("No genus names found in the provided file.")
            sys.exit(1)

        retrieve_genomes(genera_list_txt)

    else:
        # Split the input by space
        genera_space = input_argument
        genera_list_space = genera_space.split()

        if not genera_list_space:
            print("No genus names found in the provided space-separated input.")
            sys.exit(1)

        retrieve_genomes(genera_list_space)



