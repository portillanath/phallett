#!/bin/bash

# Default argument values 
script_dir=$(dirname "$0")
data_default=$(cat "$script_dir/test_genus.txt")
module=""
kmersmash=(7 9 11 12 13)
genus=""
kmersani=(12 11 10 9 8)
frag_lengths=(500)  
kmersy=(15 17 20 21 24)
kmersx=(12 11 10 9 8)
my="mash"
mx="ani"
blastpor=0.75
evalue=1e-5
file="$script_dir/GCF_000836945.fasta"
updatedb=false
all_hosts_flag=""

# Activate conda environment
conda activate enviroments

# Parse arguments
while getopts "d:m:g:a:s:f:y:x:M:X:F:u:b:e:H" opt; do
  case $opt in
    d) data_default="$OPTARG" ;;
    m) module="$OPTARG" ;;
    g) genus="$OPTARG" ;;
    a) kmersani=($OPTARG) ;;   # ANI k-mers
    s) kmersmash=($OPTARG) ;;  # Mash k-mers
    f) frag_lengths=($OPTARG) ;;
    y) kmersy=($OPTARG) ;;
    x) kmersx=($OPTARG) ;;
    M) my="$OPTARG" ;;         # Metric y-axis (mash)
    X) mx="$OPTARG" ;;         # Metric x-axis (ani)
    F) file="$OPTARG" ;;
    u) updatedb="$OPTARG" ;;
    b) blastpor="$OPTARG" ;;
    e) evalue="$OPTARG" ;;
    H) all_hosts_flag="--all-hosts" ;;
    \?) echo "Invalid option: -$OPTARG" >&2 ;;
  esac
done

# Run phallett steps as default if no module is specified
if [ -z "$module" ]; then
  bash "$script_dir/src/00.ICTV_Metadata_Resource_Resource.sh"
  bash "$script_dir/src/01A.Taxa_Curation_Level.sh" "$data_default" $all_hosts_flag
  bash "$script_dir/src/02.Bargenome.sh"
  bash "$script_dir/src/03.ANI_Metrics.sh" "${kmersmash[@]}" "$genus"
  bash "$script_dir/src/04.Mash_Metrics.sh" "${kmersani[@]}" "$genus" "${frag_lengths[@]}"
  bash "$script_dir/src/05.wraggling.sh" "${kmersx[@]}" "${kmersy[@]}" "$my" "$mx"
  bash "$script_dir/src/06.Graphing.sh" "${kmersx[@]}" "${kmersy[@]}" "$my" "$mx"
else
  case "$module" in
    ictv)
      bash "$script_dir/src/00.ICTV_Metadata_Resource.sh"
      ;;
    taxa)
      bash "$script_dir/src/01A.Taxa_Curation_Level.sh" "$data_default" $all_hosts_flag
      ;;
    file)
      bash "$script_dir/src/01B.Selecting_file.sh" "$file" "$blastpor" "$evalue" "$updatedb"
      ;;
    bargenome)
      bash "$script_dir/src/02.Bargenome.sh"
      ;;
    ani)
      bash "$script_dir/src/03.ANI_Metrics.sh" "${kmersmash[@]}" "$genus"
      ;;
    mash)
      bash "$script_dir/src/04.Mash_Metrics.sh" "${kmersani[@]}" "$genus" "${frag_lengths[@]}"
      ;;
    wraggling)
      bash "$script_dir/src/05.wraggling.sh" "${kmersx[@]}" "${kmersy[@]}" "$my" "$mx"
      ;;
    graphs)
      bash "$script_dir/src/06.Graphing.sh" "${kmersx[@]}" "${kmersy[@]}" "$my" "$mx"
      ;;
    boxplot)
      bash "$script_dir/src/07.Summary_Feed.sh"
      ;;
    *)
      echo "Unknown module: $module"
      ;;
  esac
fi

# Deactivate conda environment
conda deactivate

