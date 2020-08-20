#!/bin/bash

samples=${1:-50}
timepoints=${2:-"50 100 200 300 500 1000"}
outname=${3:-"hybrid_population"}

for t in $timepoints; do
	echo -e "$((t+1))\t1\t${samples}\t0\t${outname}.${t}.txt"
done

