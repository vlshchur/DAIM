#!/bin/bash

arrays=$1
rounds=$2
cputype=${3:-"128x24"}

if [[ -z "$arrays" || -z "$rounds" ]] ; then
    echo Must specify number of array jobs and rounds per job
    echo
    echo "Usage: $0 array_jobs rounds [cputype: 128x24]"
    exit 1
fi

sbatch -p $cputype --array=1-${arrays} run_batch.sh $rounds

echo auto_arrays=$arrays >> auto_parameters.txt
echo auto_rounds=$rounds >> auto_parameters.txt
