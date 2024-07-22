
import os
import re
import sys
import pandas as pd
from glob import glob
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import shutil
from pathlib import Path  

# Set paths
current_dir = Path.cwd() 
parent_dir = current_dir.parent
source = os.path.expanduser((Path(parent_dir) / "phallett" / "data" / "Taxa_Selected"))
subdirectories = [d for d in os.listdir(source)]
workdir = os.path.expanduser((Path(parent_dir) / "phallett" / "test" / "Metrics_Results"))

# Create output directories
for genus in subdirectories:
    genus_name = os.path.basename(os.path.normpath(genus))
    genus_dir = os.path.join(workdir, genus_name)
     # Check if the directory already exists
    if not os.path.exists(genus_dir):
       os.makedirs(genus_dir)
    #Move all the files to the corresponding directory
    if os.path.exists(genus_dir):
        files_to_move = [f for f in os.listdir(workdir) if os.path.isfile(os.path.join(workdir, f)) and genus_name in f]
        dirs_to_move= [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d)) and f"signatures_{genus_name}" in d]
       
        for file in files_to_move:
           source_path = os.path.join(workdir, file)
           destination_path = os.path.join(genus_dir, file)
    
           if os.path.exists(source_path):
              os.rename(source_path, destination_path)

        for dir in dirs_to_move:
            source_path = os.path.join(workdir, dir)
            destination_path = os.path.join(genus_dir, dir)
            if not os.path.exists(destination_path):
                os.rename(source_path, destination_path)
            else:
                print(f"Directory {destination_path} already exists")


import argparse

# Crea el analizador de argumentos
parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')

# Agrega los argumentos
parser.add_argument('-mx', type=str, help='El argumento mx')
parser.add_argument('-kmersx', type=str, help='El argumento kmersx')
parser.add_argument('-my', type=str, help='El argumento my')
parser.add_argument('-kmersy', type=str, help='El argumento kmersy')

# Analiza los argumentos
args = parser.parse_args()

# Ahora puedes acceder a los argumentos con args.mx, args.kmersx, args.my, args.kmersy
print(f'mx: {args.mx}, kmersx: {args.kmersx}, my: {args.my}, kmersy: {args.kmersy}')

mx=args.mx
my=args.my
kmersx=[int(k) for k in args.kmersx.split(",")]
kmersy=[int(k) for k in args.kmersy.split(",")]

# Set up directories on workdir
subdirectories_results = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d))]
# SEARCH FOR OUTPUT FILES
# Recover possible algorithms for each metric of the pairwise correlation
# Cases for metrics selected available
if mx == "mash":
    tool_mx = ["mash", "sourmash"]
elif mx == "ani":
    tool_mx = ["fastani", "skani"]
elif mx == "aai":
    tool_mx = ["comparem"]
elif mx == "viridic":
    tool_mx = ["viridic"]
elif mx == "vcontact2":
    tool_mx = ["vcontact2"]
else:
    raise ValueError("Invalid metric selected")

# Cases for metrics selected on Y
if my == "mash":
    tool_my = ["mash", "sourmash"]
elif my == "ani":
    tool_my = ["fastani", "skani"]
elif my == "aai":
    tool_my = ["comparem"]
elif my == "viridic":
    tool_my = ["viridic"]
elif my == "vcontact2":
    tool_my = ["vcontact2"]
else:
    raise ValueError("Invalid metric selected")

metrics = [mx, my]

# This makes the parsing for the corresponding metric
fastani_results = pd.DataFrame()
skani_results = pd.DataFrame()
mash_results =pd.DataFrame()
sourmash_results=pd.DataFrame()

for subdir in subdirectories_results:
    os.chdir(workdir)
    if os.path.exists(subdir):
       subdir_name = os.path.basename(subdir)
       if "signatures" in subdir_name:
          continue
       
       for m in metrics:
           if m == mx:
              tool_list = tool_mx
              kmers = kmersx
           if m == my:
              tool_list = tool_my
              kmers = kmersy

           for tool in tool_list:
# In case of having fastani metrics
            if tool == "fastani":
                files_fastani = glob(os.path.join(workdir, subdir_name, "fastani*"))
                files_fastani = [file for file in files_fastani if not file.endswith('.csv')]
                for file in files_fastani:
                    k = int(file.split("_")[-1].split(".")[0])
                    data_fastani = pd.read_csv(file, header=None, sep='\t')
                    data_fastani = data_fastani.iloc[:, :3]
                    data_fastani['kmer_ani'] = k
                    data_fastani['algorithm_ani'] ="fastani"
                    data_fastani.columns = ["GenomeA", "GenomeB", "ani_distance", "kmer_ani","algorithm_ani"]
                    data_fastani['GenomeA'] = data_fastani['GenomeA'].replace('.*/', '', regex=True)
                    data_fastani['GenomeB'] = data_fastani['GenomeB'].replace('.*/', '', regex=True)
                    data_fastani['GenomeA'] = data_fastani['GenomeA'].replace('.fasta', '', regex=True)
                    data_fastani['GenomeB'] = data_fastani['GenomeB'].replace('.fasta', '', regex=True)
                    fastani_results = pd.concat([fastani_results, data_fastani])
                    
            fastani_results.to_csv(os.path.join(workdir, subdir_name, f"fastani_results_{subdir_name}.csv"), index=False)
            
            # In case of having skani metrics
            if tool=="skani":
                files_skani = glob(os.path.join(workdir, subdir_name, "skani*"))
                files_skani = [file for file in files_skani if not file.endswith('.csv')]
                for file in files_skani:  
                    data_skani = pd.read_csv(file, header=None, sep='\t', names=["Ref_file","Query_file","ANI","Align_fraction_ref","Align_fraction_query","Ref_name","Query_name"])
                    data_skani = data_skani.iloc[:, :3]
                    data_skani['kmer_ani'] = "static"
                    data_skani['algorithm_ani'] ="skani"
                    data_skani = data_skani.iloc[1:]
                    data_skani.columns = ["GenomeA", "GenomeB", "ani_distance", "kmer_ani","algorithm_ani"]
                    data_skani['GenomeA'] = data_skani['GenomeA'].replace('.*/', '', regex=True)
                    data_skani['GenomeB'] = data_skani['GenomeB'].replace('.*/', '', regex=True)
                    data_skani['GenomeA'] = data_skani['GenomeA'].replace('.fasta', '', regex=True)
                    data_skani['GenomeB'] = data_skani['GenomeB'].replace('.fasta', '', regex=True)
                    skani_results = pd.concat([skani_results, data_skani])
                
            skani_results.to_csv(os.path.join(workdir, subdir_name, f"skani_results_{subdir_name}.csv"), index=False)   
            
            # In case of having mash metrics
            if tool=="mash":
                files_mash = glob(os.path.join(workdir, subdir_name, "mash*.tab"))
                for file in files_mash:
                         k = int(re.search(r'k(\d+)', file).group(1))
                         if k in kmers:
                            data_mash = pd.read_csv(file, header=None, names=["GenomeA", "GenomeB", "mash_distance", "p-value", "shared_hashes"], sep='\t')
                            data_mash = data_mash.iloc[:, :3]
                            data_mash['kmer_mash'] = k
                            data_mash['algorithm_mash'] = "mash"
                            data_mash['GenomeA'] = data_mash['GenomeA'].replace('.*/', '', regex=True)
                            data_mash['GenomeB'] = data_mash['GenomeB'].replace('.*/', '', regex=True)
                            data_mash['GenomeA'] = data_mash['GenomeA'].replace('.fasta', '', regex=True)
                            data_mash['GenomeB'] = data_mash['GenomeB'].replace('.fasta', '', regex=True)
                            mash_results = pd.concat([mash_results, data_mash])
            mash_results.to_csv(os.path.join(workdir, subdir_name, f"mash_results_{subdir_name}.csv"), index=False)   

            # In case of having sourmash metrics
            if tool=="sourmash":
              files_sourmash = glob(os.path.join(workdir, subdir_name, "sourmash*.csv"))
              accessions_path=os.path.join(source,subdir_name)
              accessions = os.listdir(accessions_path)
              accessions= [file for file in accessions if file.endswith('.fasta')]
              for file in files_sourmash:
                k_values = [int(match.group(1)) for match in re.finditer(r'k(\d+)', file)]
                k = k_values[0] if k_values else None 
                if k in kmers:
                    data_sourmash=pd.read_csv(file, sep=',') 
                    genomes= [col.split('.')[0] for col in data_sourmash.columns]  
                    data_sourmash.columns=genomes
                    data_sourmash.index=genomes
                    data_sourmash=data_sourmash.to_numpy()
                    distances=pdist(data_sourmash)
                    square_distances=squareform(distances)
                    i,j=np.triu_indices(square_distances.shape[0],k=1)
                    data_sourmash = pd.DataFrame({"GenomeA": [genomes[int(idx)] for idx in i], "GenomeB": [genomes[int(idx)] for idx in j], "mash_distance": square_distances[i.astype(int), j.astype(int)]})
                    data_sourmash['kmer_mash'] = k
                    data_sourmash['algorithm_mash']= "sourmash"
                    sourmash_results = pd.concat([sourmash_results, data_sourmash], ignore_index=True)

            sourmash_results.to_csv(os.path.join(workdir, subdir_name, f"sourmash_results_{subdir_name}.csv"), index=False)         

#Merge of differents compend
    ani_metrics_result=pd.concat([fastani_results,skani_results])
    ani_metrics_result['ani_distance'] = ani_metrics_result['ani_distance'].round(6)
    ani_metrics_result.to_csv(os.path.join(workdir,subdir_name,f"ani_metrics_{subdir_name}.csv"), index=False)
    mash_metrics_result=pd.concat([mash_results,sourmash_results])
    mash_metrics_result['mash_distance'] = mash_metrics_result['mash_distance'].round(6)
    mash_metrics_result.to_csv(os.path.join(workdir,subdir_name,f"mash_metrics_{subdir_name}.csv"), index=False)

#Create a summary of the metrics to be read as the paiwise relationships
    ani_metrics_result = ani_metrics_result.reset_index(drop=True)
    mash_metrics_result = mash_metrics_result.reset_index(drop=True)
    metrics_summary = pd.concat([ani_metrics_result, mash_metrics_result], axis=1)
    metrics_summary.rename(columns={metrics_summary.columns[4]: "algorithm_ani"}, inplace=True)
    metrics_summary.to_csv(os.path.join(workdir,subdir_name,f"summary_{subdir_name}.csv"), index=False)
    print(f"Metrics summary for {subdir_name} has been created") 
