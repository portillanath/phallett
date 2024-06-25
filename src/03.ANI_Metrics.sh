#!/bin/bash

# Default kmer values
kmers=(12 11 10 9 8)
frag_lengths=(500)  # Default fragment length
genus=""

while getopts "k:f:g:" opt; do
  case $opt in
    k)
      kmers=($OPTARG)
      ;;
    f)
      frag_lengths=($OPTARG)
      ;; 
    g)
      genus=($OPTARG)
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# Create variables for paths
source=~/phallett/data/Taxa_Selected
outdir=~/phallett/test/Metrics_Results

cd "$source"

subdirs=($(find "$source" -mindepth 1 -type d))
for subdir in "${subdirs[@]}" ; do
cd $subdir
genusname=${subdir##*/}
mv "${subdir}/signatures_${genusname}" "${outdir}/"
echo "The moving signatures were sucessful"
done

# Create a list of subdirectories within the working directory
if [ -n "$genus" ]; then
  # If a genus name is provided, only process that genus
  subdirs=("$source/$genus")
else
  # Otherwise, process all subdirectories
  subdirs=($(find "$source" -mindepth 1 -type d))
fi

# Loop through each subdirectory and create a concatenated folder
for subdir in "${subdirs[@]}"; do
  genusname=${subdir##*/}
  ls -d "${subdir}"/*.fasta > "${genusname}.list"
  mv *.list $outdir

  # Now we are going to calculate the ANI value. The file generated is interpreted as:
  # first column Ref_file, Query_file, ANI, Align_fraction_ref, Align_fraction_query, Ref_name, Query_Name

  # This calculates the triangle
  skani triangle "$source/${genusname}/"*.fasta -E > "skani_distance_${genusname}.txt"
  mv *.txt $outdir

for frag_len in "${frag_lengths[@]}"; do
    for k in "${kmers[@]}"; do
       #average_nucleotide_identity.py -i "${subdir_basename}.list" -o "${subdir}/fastani_${subdir_basename}_frag_${frag_len}_${k}" --method ANIb
      fastANI --ql "$outdir/${genusname}.list" --rl "$outdir/${genusname}.list" -o "${outdir}/fastani_${genusname}_frag_${frag_len}_${k}" --fragLen "${frag_len}" --kmer "${k}"
      if [ -f "${outdir}/fastani_${genusname}_frag_${frag_len}_${k}" ]; then
        echo "Fastani was calculated for $subdir_basename with kmer $k and fragment length $frag_len"
      else
        echo "Error: fastANI calculation failed for $subdir_basename with kmer $k and fragment length $frag_len"
      fi
    done
  done
done

