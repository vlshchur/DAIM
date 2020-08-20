batch=$1
rounds=${2:-1000}
startrounds=${3:-1}

offset=$((batch*rounds))

mkdir simulation_output/${batch}
cp input.demo.txt input.selection.txt input.output.txt simulation_output/${batch}
cd simulation_output/${batch}

subsets=$(((rounds-1)/1000))
substart=$(((startrounds-1)/1000))
for i in $(seq $substart $subsets)
do
	ii=$((i*1000))
	mkdir -p $ii
	cp input.demo.txt input.selection.txt input.output.txt $ii
done

for i in $(seq $startrounds $rounds)
do
	subbatch=$(((i-1)-(i-1)%1000))
	cd $subbatch
	sed "s/txt/${i}.txt/" input.output.txt > input.output.${i}.txt
	SELAM -d input.demo.txt -s input.selection.txt -o input.output.${i}.txt -c 1 1 -h --seed $((i + offset))
	rm input.output.${i}.txt
	cd ..
done
