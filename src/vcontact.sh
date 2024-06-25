#!/bin/bash

output_directory=~/phallett/data/Taxa_Selected
echo "Output directory: $output_directory"
source_directory=~/phallett/test/Metrics_Results
echo "Source directory: $source_directory"
vcontact_programme=~/miniconda3/envs/enviroments/bin
echo "Vcontact programme: $vcontact_programme"
# This loop generates a concatenated fasta for vcontact
subdirs=$(find "$source_directory" -type d)

# Loop through each subdirectory

for subdir in $subdirs; do
    subdir_basename=$(basename "$subdir")
    cd "$subdir"
    prodigalfile="${subdir_basename}.faa"
    prodigalpathfile="${subdir}/$prodigalfile"
    gene2genomefile="${subdir_basename}.csv"
    gene2genomepath="${subdir}/$gene2genomefile"

    cd "$vcontact_programme"

    # Run vcontact2 gene2genome
    ./vcontact2_gene2genome -p "$prodigalpathfile" -o "${subdir}/$gene2genomefile" -s 'Prodigal-FAA'

    # Running the complete workflow
    ./vcontact2 --raw-proteins "$prodigalpathfile" --rel-mode 'Diamond' --proteins-fp "$gene2genomepath" --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /hpcfs/home/ciencias_biologicas/na.portilla10/.conda/envs/vcontact2/bin/cluster_one-1.0.jar --output-dir "$subdir" -t 8
done
