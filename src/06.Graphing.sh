#!/bin/bash
#Author:Nathalia Portilla
#name script:wraggling.sh
#SBATCH -o Graphing_job.o%j

kmersy=(7,9,11,12,13)
kmersx=(12,11,10,9,8)
my="mash"
mx="ani"
parent_dir=$(dirname "$PWD")

while getopts "mx:kmersx:my:kmersy" option; do
    case $option in
        mx) #Handle the -mx flag with an argument
        mx=${OPTARG}
        ;;
        kmersx) #Handle the -kmersx flag with an argument
        kmersx=${OPTARG} 
        ;;
        my) #Handle the -my flag with an argument
        my=${OPTARG}
        ;;
        kmersy) #Handle the -kmersy flag with an argumentcd
        kmersy=${OPTARG}
        ;;
    esac
done

python3 "$parent_dir/phallett/src/06.Graphing.py" -mx "$mx" -kmersx "$kmersx" -my "$my" -kmersy "$kmersy" > Graphing.log 2>&1
cat Graphing.log
