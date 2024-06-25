#!/bin/bash
#Default argument values 
data_default=~/phallett/test_genus.txt
module=""
kmersmash=(7 9 11 12 13)
genus=""
kmersani=(12 11 10 9 8)
frag_lengths=(500)  
kmersy=(15,17,20,21,24)
kmersx=(12,11,10,9,8)
my="mash"
mx="ani"
blastpor=0.75
evalue=1e-5
file="~/phallett/GCF_000836945.fasta"
updatedb=false

conda activate enviroments 

#Parse arguments
while getopts "d:m:g:ka:km:f:ky:kx:my:mx:fl:u:b:e" opt; do
  case $opt in
    d)
      data_default=$OPTARG
      ;;
    m)
      module=$OPTARG
      ;;
    g)
      genus=$OPTARG
      ;;
    ka)
      kmersani=$OPTARG
      ;;
    km)
      kmersmash=$OPTARG
      ;;
    f)
      frag_lengths=$OPTARG
      ;;
    ky)
      kmersy=$OPTARG
      ;;
    kx)
      kmersx=$OPTARG
      ;;
    my)
      my=$OPTARG
      ;;
    mx)
      mx=$OPTARG
      ;;
    fl)
      file=$OPTARG
      ;;
    u)  
      updatedb=$OPTARG
      ;;
    b)
      blastpor=$OPTARG
      ;;
    e)
      evalue=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

#Run phallet steps with flags
if [ -z "$module" ]; then
  bash ~/phallett/src/00.ICTV_Metadata_Resource.sh
  bash ~/phallett/src/01.Taxa_Curation_Level.sh $data_default
  bash ~/phallett/src/02.Mash_Metrics.sh $kmersmash $genus
  bash ~/phallett/src/03.ANI_Metrics.sh $kmersani $genus $frag_lengths
  bash ~/phallett/src/04.wraggling.sh $kmersx $kmersy $my $mx
  bash ~/phallett/src/05.Graphing.sh $kmersx $kmersy $my $mx
else
  if [ "$module" == "ictv" ]; then
    bash ~/phallett/src/00.ICTV_Metadata_Resource.sh
  elif [ "$module" == "taxa" ]; then
    bash ~/phallett/src/01.Taxa_Curation_Level.sh $data_default
  elif [ "$module" == "file" ]; then
    bash ~/phallett/src/01B.Selecting_file.sh $file $blastpor $evalue $updatedb
  elif [ "$module" == "mash" ]; then
    bash ~/phallett/src/02.Mash_Metrics.sh $kmersmash $genus
  elif [ "$module" == "ani" ]; then
    bash ~/phallett/src/03.ANI_Metrics.sh $kmersani $genus $frag_lengths
  elif [ "$module" == "wraggling" ]; then
    bash ~/phallett/src/04.wraggling.sh $kmersx $kmersy $my $mx
  elif [ "$module" == "graphs" ]; then
    bash ~/phallett/src/05.Graphing.sh $kmersx $kmersy $my $mx
  fi
fi

#Deactivate conda enviroments
conda deactivate
