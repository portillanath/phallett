import argparse
import os
import re
import subprocess
import glob
import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import urllib.error

#Runninf selecting file
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-file', type=str)
parser.add_argument('-updatedb', type=str, default='False')
parser.add_argument('-blastpor', type=float)
parser.add_argument('-evalue', type=float)
args = parser.parse_args()
#Convert args.updatedb to boolean
updatedb = False if args.updatedb.lower() == 'false' else True
print(vars(args))

if updatedb:
    ICTV_database = os.path.expanduser("~/phallett/data/ICTV_database")
    os.makedirs(ICTV_database, exist_ok=True)
    Entrez.email = "na.portilla10@uniandes.edu.co"
    guide_accessions= pd.read_csv(os.path.expanduser("~/phallett/data/Virus_Metadata_Resource/VMR.csv"))

    for accession in guide_accessions['Virus GENBANK accession']:
        try:
            print(f"Downloading {accession}")
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            with open(os.path.join(ICTV_database, f"{accession}.fasta"), "w") as f:
                f.write(handle.read())
        except Exception as e:
            print(f"Error al descargar {accession}: {e}. Saltando al siguiente.")
            continue
else:
    print("Database not updated")

if args.file:
    #Create a new analysis folder 
    existing_folders=glob.glob(os.path.expanduser("~/phallett/data/Taxa_Selected/*"))
    existing_numbers = [int(re.search(r'analysis_(\d+)', folder).group(1)) for folder in existing_folders if re.search(r'analysis_(\d+)', folder)]
    next_number=max(existing_numbers)+1 if existing_numbers else 1
    analysis_folder=os.path.expanduser(f"~/phallett/data/Taxa_Selected/analysis_{next_number}")
    os.makedirs(analysis_folder, exist_ok=True)
    print(f"Analysis folder created: {analysis_folder}")
    #Read the fasta file
    file_path = os.path.expanduser('~/phallett/GCF_000836945.fasta')
    fasta_string = open(file_path).read()
    print(fasta_string)
    
    #Format the database
    command = "cat ~/phallett/data/ICTV_database/*.fasta > ~/phallett/data/ICTV_database_combined.fasta"
    subprocess.run(command, shell=True)
    command="makeblastdb -in ~/phallett/data/ICTV_database_combined.fasta -dbtype nucl" 
    makeblastdb=subprocess.run(command, shell=True)  

    #Save the query fasta file
    fasta_file_name = os.path.basename(args.file)  # Get the original file name
    fasta_file_path = os.path.join(analysis_folder, fasta_file_name)
    with open(fasta_file_path, "w") as f:
        f.write(fasta_string)
    print("Query Fasta file saved")

    #Search on Blast
    command = f"blastn -query {fasta_file_path} -db ~/phallett/data/ICTV_database_combined.fasta -out result.xml -outfmt 5"  # Apunta a la base de datos combinada
    result= subprocess.run(command, shell=True)
    
    if result.returncode == 0:
        print("BLAST executed successfully")
    else:
        print("Error al ejecutar BLAST")
 
# Parse BLAST results
result_handle = open("result.xml")
blast_records = NCBIXML.parse(result_handle)

def clean_filename(filename):
    # Mantén todo el texto después del primer espacio
    filename = filename.split(' ', 1)[-1]
    # Borra todo después de un punto
    filename = filename.split('.')[0]
    # Elimina los caracteres no válidos
    invalid_chars = ['<', '>', ':', '"', '/', '\\', '|', '?', '*']
    for char in invalid_chars:
        filename = filename.replace(char, '')
    return filename

#Dowload sequences from blast
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        print(f"Alignment: {alignment.title}")
        for hsp in alignment.hsps:
            #Filter by evalue and identity
            if hsp.identities/len(hsp.query) >= args.blastpor and hsp.expect < args.evalue:
                print(f"Found match: {alignment.title}")
                from Bio import Entrez
            if alignment.accession:
                try:
                   handle = Entrez.efetch(db="nucleotide", id=alignment.accession, rettype="fasta", retmode="text")
                except urllib.error.HTTPError as err:
                  print(f"Error doing efetch: {err}")
            else:
                 print("Accession ID is empty. Cannot do efetch request.")
                
            # Clean the alignment title to make it a valid filename
            valid_filename = clean_filename(alignment.title)
            # Now you can open the file with the new name
            with open(os.path.join(analysis_folder, f"{valid_filename}.fasta"), "w") as f:
                f.write(handle.read())
#remove empty files
def remove_empty_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            file_path = os.path.join(directory, filename)
            if os.path.getsize(file_path) == 0:
                print(f"Removing empty file: {file_path}")
                os.remove(file_path)

# Call the function on the analysis folder
remove_empty_files(analysis_folder)