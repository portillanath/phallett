#!/bin/bash

blastpor=0.75
evalue=1e-5
parent_dir=$(dirname "$PWD")
file="$parent_dir/phallett/GCF_000836945.fasta"
updatedb="true"

while getopts "fl:bp:e:u" opt; do
  case $opt in
    f)
      file=$OPTARG
      ;;
    b)
      blastpor=$OPTARG
      ;;
    e)
      evalue=$OPTARG
      ;;
    u)
      updatedb=$OPTARG
      ;;
  esac
done

#Now is running 
python3 $parent_dir/phallett/src/01B.Selecting_file.py -file $file -updatedb $updatedb -blastpor $blastpor -evalue $evalue
 
