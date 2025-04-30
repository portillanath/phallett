#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import pandas as pd
from Bio import Entrez
from pathlib import Path

def retrieve_genomes(genera_list, all_hosts=False):
    # Setup paths
    current_dir = Path.cwd()
    parent_dir = current_dir.parent
    print(f"Working in: {parent_dir}")

    # Load and clean VMR data
    vmr_path = parent_dir / "phallett/data/Virus_Metadata_Resource/VMR.csv"
    ICTV_assignation = pd.read_csv(vmr_path)
    ICTV_assignation = ICTV_assignation.rename(columns=lambda x: x.strip())
    ICTV_assignation.columns = [
        "Sort","Isolate_Sort","Realm","Subrealm","Kingdom","Subkingdom","Phylum","Subphylum","Class","Subclass",
        "Order","Suborder","Family","Subfamily","Genus","Subgenus","Species","Exemplar_or_additional_isolate",
        "Virus_name(s)","Virus_name_abbreviation(s)","Virus_isolate_designation","Virus_GENBANK_accession",
        "Virus_REFSEQ_accession","Genome_coverage","Genome_composition","Host_source"
    ]

    # Filter based on coverage
    ICTV_assignation = ICTV_assignation[
        ICTV_assignation["Genome_coverage"].isin(["Complete genome", "Complete coding genome"])
    ]

    # Optional filter by host source
    if not all_hosts:
        ICTV_assignation = ICTV_assignation[
            ICTV_assignation["Host_source"].isin(["archaea", "bacteria"])
        ]

    # Output filtered metadata
    filtered_csv = parent_dir / "phallett/test/ICTV_assignation_Complete.csv"
    ICTV_assignation.to_csv(filtered_csv, index=False)

    # Get unique accessions
    accessions = ICTV_assignation["Virus_GENBANK_accession"].dropna().unique()
    print(f"Number of phages with complete genome: {len(accessions)}")

    # Create output folder
    ncbi_genome_actual = os.path.expanduser(str(parent_dir / "phallett/data/Taxa_Selected"))
    os.makedirs(ncbi_genome_actual, exist_ok=True)

    for genus in genera_list:
        genus_data = ICTV_assignation[ICTV_assignation["Genus"] == genus]

        if genus_data.empty:
            print(f"No records found for the genus: {genus}")
            continue

        genus_folder = os.path.join(ncbi_genome_actual, genus)
        os.makedirs(genus_folder, exist_ok=True)

        for accession in genus_data["Virus_GENBANK_accession"]:
            if pd.isna(accession):
                continue
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
                    print(f"Error fetching {accession}: {e}")
                    time.sleep(5)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: script.py <genus_file.txt or space-separated genus> [--all-hosts]")
        sys.exit(1)

    input_argument = sys.argv[1]
    all_hosts_flag = "--all-hosts" in sys.argv

    if os.path.exists(input_argument):
        with open(input_argument, 'r') as file:
            genera_list_txt = [line.strip() for line in file if line.strip()]
        if not genera_list_txt:
            print("No genus names found in the provided file.")
            sys.exit(1)
        retrieve_genomes(genera_list_txt, all_hosts=all_hosts_flag)
    else:
        genera_list_space = input_argument.split()
        if not genera_list_space:
            print("No genus names found in space-separated input.")
            sys.exit(1)
        retrieve_genomes(genera_list_space, all_hosts=all_hosts_flag)



