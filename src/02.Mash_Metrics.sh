#!/bin/bash

#Default kmers values
kmers=(7 9 11 12 13)
genus=""

while getopts "k:g:" opt;do
    case ${opt} in
        k )
            kmers=(${OPTARG})
            ;;
        g ) 
            genus=(${OPTARG})
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            ;;
    esac
done

# Create variables for paths
parent_dir=$(dirname "$PWD")
source="$parent_dir/phallett/data/Taxa_Selected"
mkdir -p"$parent_dir/phallett/test/Metrics_Results"
outdir="$parent_dir/phallett/test/Metrics_Results"

# Create a list of subdirectories within the working directory
if [ -n "$genus" ]; then
  # If a genus name is provided, only process that genus
  subdirs=("$source/$genus")
else
  # Otherwise, process all subdirectories
  subdirs=($(find "$source" -mindepth 1 -type d))
fi

# Create a single signatures directory
sig_subdir_basename="signatures"
mkdir -p "${outdir}/${sig_subdir_basename}"

# Loop through each subdirectory
for subdir in "${subdirs[@]}"; do
  subdir_basename=$(basename "$subdir")
  cd "$subdir"

  # MASH ALGORITHM
  # Create a sketch of all the sequences
  # Use 64-bit hashes and a sketch default of 1000
  for k in "${kmers[@]}"; do
    sketch_file="${subdir}/sketch_mash_${subdir_basename}_${k}"
    if ls ${subdir}/*.fasta 1> /dev/null 2>&1; then
    mash sketch -k ${k} -o ${sketch_file} ${subdir}/*.fasta
    echo "The sketching is complete for $subdir_basename with kmer $k"
    # Calculate mash distance with the sketches
    if [ -f "${sketch_file}.msh" ]; then
      mash dist ${sketch_file}.msh ${sketch_file}.msh > ${subdir}/mash_distance_${subdir_basename}_k${k}.tab
      echo "The matrix distance is calculated for $subdir_basename with kmer $k"
    else
      echo "Sketch file ${sketch_file}.msh does not exist"
    fi
    else
      echo "No fasta files found in ${subdir}"
    fi
  done

  mv sketch* $outdir
  mv *.tab $outdir

  #SOURMASH ALGORITHM

  # Create an "," array of kmers 
  kmerslist=""
  for k in "${kmers[@]}";do
    kmerslist="${kmerslist}${k},"
  done
  kmerslist=${kmerslist::-1}
  echo "$kmerslist"
 
  sourmash compute --ksizes "${kmerslist}" *.fasta --singleton

  #Move the signatures to the newly created subdirectory
  mv *.sig "${outdir}/${sig_subdir_basename}/"
  echo "All the signatures were moved"

#Run sourmash compare inside the current subdir
for k in ${kmers[@]}; do
  sourmash compare --k ${k} --csv ${outdir}/sourmash_distance_${subdir_basename}_k${k}.csv ${outdir}/${sig_subdir_basename}/*.sig
  echo "The matrix distance is calculated for $subdir_basename with kmer $k"
done

mv sourmash* $outdir
 
done