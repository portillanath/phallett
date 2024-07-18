import argparse
import os
import re
import subprocess
import glob
import pandas as pd
from Bio.Blast import NCBIXML
from Bio import Entrez
import urllib.error
from pathlib import Path  # Corrected import statement

# Argument parsing
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-file', type=str, help='Path to the fasta file')
parser.add_argument('-updatedb', type=str, default='False', help='Whether to update the database')
parser.add_argument('-blastpor', type=float, help='BLAST percentage identity threshold')
parser.add_argument('-evalue', type=float, help='BLAST e-value threshold')
args = parser.parse_args()

# Convert `updatedb` argument to boolean
updatedb = args.updatedb.lower() == 'true'

print(vars(args))

if updatedb:
    current_dir = Path.cwd()  # Get the current working directory
    parent_dir = current_dir.parent
    ICTV_database = Path(parent_dir) / "phallett" / "data" / "ICTV_database"
    ICTV_database.mkdir(parents=True, exist_ok=True)
    Entrez.email = "na.portilla10@uniandes.edu.co"
    guide_accessions = pd.read_csv(Path(parent_dir) / "phallett" / "data" / "Virus_Metadata_Resource" / "VMR.csv")

    for accession in guide_accessions['Virus GENBANK accession']:
        try:
            print(f"Downloading {accession}")
            with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
                with open(ICTV_database / f"{accession}.fasta", "w") as f:
                    f.write(handle.read())
        except Exception as e:
            print(f"Error al descargar {accession}: {e}. Saltando al siguiente.")

else:
    print("Database not updated")

if args.file:
    # Create a new analysis folder
    existing_folders = glob.glob(str(Path(parent_dir) / "phallett" / "data" / "Taxa_Selected" / "*"))
    existing_numbers = [int(re.search(r'analysis_(\d+)', folder).group(1)) for folder in existing_folders if re.search(r'analysis_(\d+)', folder)]
    next_number = max(existing_numbers) + 1 if existing_numbers else 1
    analysis_folder = Path(parent_dir) / "phallett" / "data" / "Taxa_Selected" / f"analysis_{next_number}"
    analysis_folder.mkdir(parents=True, exist_ok=True)
    print(f"Analysis folder created: {analysis_folder}")

    # Read the fasta file
    file_path = Path(parent_dir) / "phallett" / "GCF_000836945.fasta"
    fasta_string = file_path.read_text()
    print(fasta_string)

    # Format the database
    combined_fasta = Path(parent_dir) / "phallett" / "data" / "ICTV_database_combined.fasta"
    command = f"cat {ICTV_database}/*.fasta > {combined_fasta}"
    subprocess.run(command, shell=True)

    # Create the BLAST database
    command = f"makeblastdb -in {combined_fasta} -dbtype nucl"
    subprocess.run(command, shell=True)

    # Save the query fasta file
    fasta_file_name = os.path.basename(args.file)  # Get the original file name
    fasta_file_path = analysis_folder / fasta_file_name
    fasta_file_path.write_text(fasta_string)
    print("Query Fasta file saved")

    # Search on BLAST
    command = f"blastn -query {fasta_file_path} -db {combined_fasta} -out result.xml -outfmt 5"
    result = subprocess.run(command, shell=True)

    if result.returncode == 0:
        print("BLAST executed successfully")
    else:
        print("Error al ejecutar BLAST")

# Parse BLAST results
result_handle = open("result.xml")
blast_records = NCBIXML.parse(result_handle)

def clean_filename(filename):
    # Clean the filename
    filename = filename.split(' ', 1)[-1]
    filename = filename.split('.')[0]
    invalid_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*']
    for char in invalid_chars:
        filename = filename.replace(char, '')
    return filename

# Download sequences from BLAST
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        print(f"Alignment: {alignment.title}")
        for hsp in alignment.hsps:
            if hsp.identities / len(hsp.query) >= args.blastpor and hsp.expect < args.evalue:
                print(f"Found match: {alignment.title}")
                if alignment.accession:
                    try:
                        with Entrez.efetch(db="nucleotide", id=alignment.accession, rettype="fasta", retmode="text") as handle:
                            valid_filename = clean_filename(alignment.title)
                            with open(analysis_folder / f"{valid_filename}.fasta", "w") as f:
                                f.write(handle.read())
                    except urllib.error.HTTPError as err:
                        print(f"Error doing efetch: {err}")
                else:
                    print("Accession ID is empty. Cannot do efetch request.")

# Remove empty files
def remove_empty_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            file_path = os.path.join(directory, filename)
            if os.path.getsize(file_path) == 0:
                print(f"Removing empty file: {file_path}")
                os.remove(file_path)

# Call the function on the analysis folder
remove_empty_files(analysis_folder)
