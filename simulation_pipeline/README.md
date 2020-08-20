Brief instructions for running simulations as used in Shchur et al (2020).

Install SELAM from https://github.com/russcd/SELAM and compiled as described in the Github repository.
Place the SELAM binary in your $PATH.
Copy all files and subdirectories in this directory to you working directory.
Run `make_sim_dir.sh` to generate all SELAM input files and place them in a new directory. Example:
```
$  ./make_sim_dir.sh 0.1 0.01 10000 200 "100 200 1000"
```

This will set up the simulations to run a scenario with a an introgression pulse of 0.1 (10%), a diploid homozygous selection advantage of 0.01, with a population size of 10000 diploid individuals. 200 individuals will be sampled from each simulation at 100, 200 and 1000 generations. A directory named `admix_0.1-sel_0.01/` will be created. If it allready exists, you will get an error message.

Enter this director. To start a single batch of simulations, run the `batch_run.sh` script. Example:

```
./batch_sim.sh 0 1000 1
```

This will start a batch of 1000 simulations. The batch has batch number "0" and will start on simulation #1  and go up to simulation #1000. This simulation number is used as a seed for the random number generation and will ensure that each of the 1000 simulations are unique, but that they can be rerun (using the same seed number if necessary). Each simulation will be run sequentially, and if the batch is aborted (after for instance 50 simulations, you can restart at simulation 50 by restarting them using the following command `./batch_sim.sh 0 1000 50`.

In order to run many more simulations in parallel, you can for instance GNU parallel and start several instances of `batch_sim.sh`. Example:

```
parallel ./batch_sim.sh {} 1000 1 ::: {1..10}
```

will start 10 instances in parallel. (Though take note that you shouldn't start more parallel jobs than you have available processor threads/cores.)

If you are running your simulations on a cluster using SLURM, an example SLURM script is provided in `run_batch.sh`, though note that it will probably have to be edited to suit the exact SLURM configuration used on your cluster. In uses the `--array` function that some SLURM setups support, and in order to for instance run 10 batches of 1000 simulations you can use the following command:

```
sbatch --array=1-10 ./run_batch.sh 1000
```

By default it will run 72 jobs.

The simulations will be saved in `simulation_output/`, and in subdirectories separated by batch ID, and then further in subdirectories containing 1000 simulations each. This is done to prevent problems that can be caused in Linux by having too many files in a single directory. For instance, the output of simulations 1000-1999 in batch 20 are found in `simulation_output/20/1000/`.

SELAM outputs haplotypes, specified by genotype identity (in a simple 2-population case 0 or 1) and the length of the haplotypes. Based on this, it is easy to calculate genotype frequencies and tract length using awk. For example, in all simulation used in our study, the selected allele comes from genotype 1 and is located at position 0.5 (the middle of the chromosome). For example, to calculate the average genotype frequency and the standard deviation of all simulations sampled at timepoint 1000 generations, the following command can be used:

```
for x in simulation_output/*/*/hybrid_population.1000.*.txt; do cat $x | awk '$1 == "##" {j+=1} $7 == 1 && $8 < 0.5 && $9 > 0.5 {i+=1} END {print i/j}'; done > genotype_frequencies.txt
```

where `$1 == "##"` counts the number of samples and `$7 == 1 && $8 < 0.5 && $9 > 0.5` identifies haplotype blocks carrying the selected site. To calculate haplotype lengths:

```
for x in simulation_output/*/*/hybrid_population.1000.*.txt; do cat $x | awk '$7 == 1 && $8 < 0.5 && $9 > 0.5 {print $9-$8; exit}'; done > haplotype_lengths.txt
```

Since the haplotype lengths in simulations are not independent, we only sample the first selected haplotype of each simulation.

