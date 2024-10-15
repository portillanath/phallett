import os
import re
import pandas as pd
from glob import glob
import numpy as np
from scipy.spatial.distance import pdist, squareform
from pathlib import Path  
import argparse

# Set paths
current_dir = Path.cwd() 
parent_dir = current_dir.parent
source = os.path.expanduser((Path(parent_dir) / "phallett" / "data" / "Taxa_Selected"))
workdir = os.path.expanduser((Path(parent_dir) / "phallett" / "test" / "Metrics_Results"))

# List subdirectories in the source directory
subdirectories = [d for d in os.listdir(source) if os.path.isdir(os.path.join(source, d))]

# Create output directories and move files
for genus in subdirectories:
    genus_name = os.path.basename(os.path.normpath(genus))
    genus_dir = os.path.join(workdir, genus_name)
    
    # Check if the directory already exists
    if not os.path.exists(genus_dir):
        os.makedirs(genus_dir)
    
    # Move files and directories
    files_to_move = [f for f in os.listdir(workdir) if os.path.isfile(os.path.join(workdir, f)) and genus_name in f]
    dirs_to_move = [d for d in os.listdir(workdir) if os.path.isdir(os.path.join(workdir, d)) and f"signatures_{genus_name}" in d]

    for file in files_to_move:
        source_path = os.path.join(workdir, file)
        destination_path = os.path.join(genus_dir, file)
        if os.path.exists(source_path):
            os.rename(source_path, destination_path)

    for dir in dirs_to_move:
        source_path = os.path.join(workdir, dir)
        destination_path = os.path.join(genus_dir, dir)
        if not os.path.exists(destination_path):
            shutil.move(source_path, destination_path)
        else:
            print(f"Directory {destination_path} already exists")

# Create the argument parser
parser = argparse.ArgumentParser(description='Process some arguments.')
parser.add_argument('-mx', type=str, help='The mx argument')
parser.add_argument('-kmersx', type=str, help='The kmersx argument')
parser.add_argument('-my', type=str, help='The my argument')
parser.add_argument('-kmersy', type=str, help='The kmersy argument')

# Parse the arguments
args = parser.parse_args()

mx = args.mx
my = args.my
kmersx = [int(k) for k in args.kmersx.split(",")]
kmersy = [int(k) for k in args.kmersy.split(",")]

# Define metrics based on input arguments
metric_tools = {
    "mash": ["mash", "sourmash"],
    "ani": ["fastani", "skani"],
    "aai": ["comparem"],
    "viridic": ["viridic"],
    "vcontact2": ["vcontact2"]
}

tool_mx = metric_tools.get(mx)
tool_my = metric_tools.get(my)

if tool_mx is None or tool_my is None:
    raise ValueError("Invalid metric selected")

metrics = [mx, my]

# Process each genus directory
for subdir in subdirectories:
    genus_name = os.path.basename(os.path.normpath(subdir))
    genus_dir = os.path.join(workdir, genus_name)

    # Reinitialize DataFrames for each genus
    fastani_results = pd.DataFrame()
    skani_results = pd.DataFrame()
    mash_results = pd.DataFrame()
    sourmash_results = pd.DataFrame()

    # Change to the genus directory
    os.chdir(genus_dir)

    for m in metrics:
        if m == mx:
            tool_list = tool_mx
            kmers = kmersx
        else:
            tool_list = tool_my
            kmers = kmersy

        for tool in tool_list:
            if tool == "fastani":
                files_fastani = glob(os.path.join(genus_dir, "fastani*"))
                files_fastani = [file for file in files_fastani if not file.endswith('.csv')]
                for file in files_fastani:
                    k = int(file.split("_")[-1].split(".")[0])
                    data_fastani = pd.read_csv(file, header=None, sep='\t')
                    data_fastani = data_fastani.iloc[:, :3]
                    data_fastani['kmer_ani'] = k
                    data_fastani['algorithm_ani'] = "fastani"
                    data_fastani.columns = ["GenomeA", "GenomeB", "ani_distance", "kmer_ani", "algorithm_ani"]
                    data_fastani['GenomeA'] = data_fastani['GenomeA'].replace('.*/', '', regex=True).replace('.fasta', '', regex=True)
                    data_fastani['GenomeB'] = data_fastani['GenomeB'].replace('.*/', '', regex=True).replace('.fasta', '', regex=True)
                    fastani_results = pd.concat([fastani_results, data_fastani], ignore_index=True)

                fastani_results.to_csv(os.path.join(genus_dir, f"fastani_results_{genus_name}.csv"), index=False)
            
            elif tool == "skani":
                files_skani = glob(os.path.join(genus_dir, "skani*"))
                files_skani = [file for file in files_skani if not file.endswith('.csv')]
                for file in files_skani:
                    data_skani = pd.read_csv(file, header=None, sep='\t', names=["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query", "Ref_name", "Query_name"])
                    data_skani = data_skani.iloc[:, :3]
                    data_skani['kmer_ani'] = "static"
                    data_skani['algorithm_ani'] = "skani"
                    data_skani = data_skani.iloc[1:]
                    data_skani.columns = ["GenomeA", "GenomeB", "ani_distance", "kmer_ani", "algorithm_ani"]
                    data_skani['GenomeA'] = data_skani['GenomeA'].replace('.*/', '', regex=True).replace('.fasta', '', regex=True)
                    data_skani['GenomeB'] = data_skani['GenomeB'].replace('.*/', '', regex=True).replace('.fasta', '', regex=True)
                    skani_results = pd.concat([skani_results, data_skani], ignore_index=True)
                
                skani_results.to_csv(os.path.join(genus_dir, f"skani_results_{genus_name}.csv"), index=False)
            
            elif tool == "mash":
                files_mash = glob(os.path.join(genus_dir, "mash*.tab"))
                for file in files_mash:
                    k = int(re.search(r'k(\d+)', file).group(1))
                    if k in kmers:
                        data_mash = pd.read_csv(file, header=None, names=["GenomeA", "GenomeB", "mash_distance", "p-value", "shared_hashes"], sep='\t')
                        data_mash = data_mash.iloc[:, :3]
                        data_mash['kmer_mash'] = k
                        data_mash['algorithm_mash'] = "mash"
                        data_mash['GenomeA'] = data_mash['GenomeA'].replace('.*/', '', regex=True).replace('.fasta', '', regex=True)
                        data_mash['GenomeB'] = data_mash['GenomeB'].replace('.*/', '', regex=True).replace('.fasta', '', regex=True)
                        mash_results = pd.concat([mash_results, data_mash], ignore_index=True)

                mash_results.to_csv(os.path.join(genus_dir, f"mash_results_{genus_name}.csv"), index=False)

            elif tool == "sourmash":
                files_sourmash = glob(os.path.join(genus_dir, "sourmash*.csv"))
                accessions_path = os.path.join(source, genus_name)
                accessions = [file for file in os.listdir(accessions_path) if file.endswith('.fasta')]
                for file in files_sourmash:
                    k_values = [int(match.group(1)) for match in re.finditer(r'k(\d+)', file)]
                    k = k_values[0] if k_values else None
                    if k in kmers:
                        data_sourmash = pd.read_csv(file, sep=',')
                        genomes = [col.split('.')[0] for col in data_sourmash.columns]
                        data_sourmash.columns = genomes
                        data_sourmash.index = genomes
                        data_sourmash = data_sourmash.to_numpy()
                        distances = pdist(data_sourmash)
                        square_distances = squareform(distances)
                        i, j = np.triu_indices(square_distances.shape[0], k=1)
                        data_sourmash = pd.DataFrame({
                            "GenomeA": [genomes[int(idx)] for idx in i],
                            "GenomeB": [genomes[int(idx)] for idx in j],
                            "mash_distance": square_distances[i.astype(int), j.astype(int)]
                        })
                        data_sourmash['kmer_mash'] = k
                        data_sourmash['algorithm_mash'] = "sourmash"
                        sourmash_results = pd.concat([sourmash_results, data_sourmash], ignore_index=True)

                sourmash_results.to_csv(os.path.join(genus_dir, f"sourmash_results_{genus_name}.csv"), index=False)

    # Merge results
    ani_metrics_result = pd.concat([fastani_results, skani_results], ignore_index=True)
    ani_metrics_result['ani_distance'] = ani_metrics_result['ani_distance'].round(6)
    ani_metrics_result.to_csv(os.path.join(genus_dir, f"ani_metrics_{genus_name}.csv"), index=False)

    mash_metrics_result = pd.concat([mash_results, sourmash_results], ignore_index=True)
    mash_metrics_result['mash_distance'] = mash_metrics_result['mash_distance'].round(6)
    mash_metrics_result.to_csv(os.path.join(genus_dir, f"mash_metrics_{genus_name}.csv"), index=False)

    # Create a summary of the metrics
    metrics_summary = pd.merge(ani_metrics_result, mash_metrics_result, on=['GenomeA', 'GenomeB'], how='inner')
    metrics_summary.to_csv(os.path.join(genus_dir, f"summary_{genus_name}.csv"), index=False)
    print(f"Metrics summary for {genus_name} has been created")
