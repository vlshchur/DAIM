#!/bin/bash

DEFDIR="./default_directory/"

admix=$1
sel=$2
pop=${3:-100000}
samples=${4:-50}
timepoints=${5:-"50 100 200 300 500 1000"}
outname=${6:-"hybrid_population"}

if [[ -z "$admix" || -z "$sel" ]] ; then
    echo "Must specify admixture fraction and selective coeffient"
    echo
    echo "Usage: $0 admixture selection [population_size, number_of_samples, timepoints, output_filename]"
    exit 1
fi

dirname=admix_${admix}-sel_${sel}

if [ $pop -ne 10000 ]
then
	dirname=${dirname}-pop_${pop}
fi

if [ -d $dirname ]
then
	echo Directory $dirname exists!
	exit 1
fi

cp -r $DEFDIR ${dirname}

if (( $(echo "$sel != 0" |bc -l) ))
then
	./generate_selection_file.py d ${sel} > ${dirname}/input.selection.txt
fi

./generate_demo_file.py ${admix} ${pop} > ${dirname}/input.demo.txt
./generate_sim_output.sh $samples "$timepoints" $outname > ${dirname}/input.output.txt

cp ${dirname}/input.* ${dirname}/simulation_output/

echo auto_admix=$admix >> ${dirname}/auto_parameters.txt
echo auto_sel=$sel >> ${dirname}/auto_parameters.txt
echo auto_pop=$pop >> ${dirname}/auto_parameters.txt
echo auto_samples=$samples >> ${dirname}/auto_parameters.txt
echo auto_timepoints=\"$timepoints\" >> ${dirname}/auto_parameters.txt
echo auto_outname=$outname >> ${dirname}/auto_parameters.txt

