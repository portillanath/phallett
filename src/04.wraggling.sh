#!/bin/bash
#Author:Nathalia Portilla
#name script:wraggling.sh
#SBATCH -o Metrics_Wraggling_job.o%j

kmersy=(7,9,11,12,13)
kmersx=(12,11,10,9,8)
my="mash"
mx="ani"

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
        kmersy) #Handle the -kmersy flag with an argument
        kmersy=${OPTARG}
        ;;
    esac
done

python3 ~/phallett/src/04.wraggling.py -mx "$mx" -kmersx "$kmersx" -my "$my" -kmersy "$kmersy" > Metrics_Wraggling.log 2>&1
cat Metrics_Wraggling.log


