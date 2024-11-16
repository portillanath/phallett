#!/bin/bash
#Author:Nathalia Portilla
#name script:summary_feed.sh
#SBATCH -o Metrics_Wraggling_job.o%j

parent_dir=$(dirname "$PWD")
python3 "$parent_dir/src/08.Summary_Feed.py"

