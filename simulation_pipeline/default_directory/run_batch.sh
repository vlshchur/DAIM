#!/bin/bash
#
#SBATCH --job-name=selection
#SBATCH --output=./sbatch_output/%a.out
#SBATCH --error=./sbatch_output/%a.err
#SBATCH --ntasks=1
#SBATCH -p 128x24
#SBATCH --mem 5G
#SBATCH --exclusive=user
#SBATCH --array=1-72


module load gnu7/7.3.0

rounds=${1:-1400}
startrounds=${2:-1}

#Do this 
./batch_sim.sh ${SLURM_ARRAY_TASK_ID} $rounds $startrounds
