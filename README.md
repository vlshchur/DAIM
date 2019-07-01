# tract_length
Calculator for the tract length distribution under Deterministic Adaptive Introgression Model (DAIM).

## Model trajectory with logistic function
Choose mode `./tract_length logit`, then add at least two points with `--point/-p` argument. Example:
`./tract_length.py logit -p 0 None 0.01 -p 1500 0.005 0`
`--points/-p` is followed by the time (measured in generations, backward, so that `0` is present time), allele frequency (AF) at this time and selection coefficient. One parameter should be set to `None`, and will be calculated from the corresponding logistic function (with an exception of neutral case, where you just need to make sure that AFs at the start and at the end are the same). Simulation with the period of neutrality:
`./tract_length.py logit -p 0 None 0.01 -p 1500 0.005 0 -p 2000 0.005 0`

## Model trajectory with precomputed function
Choose mode `./tract_length precomp` followed by at least one input file.
Input file is tab-seperated. The lines with expected (mean) trajectories should start with `MT`. Then the AF is given for each generation starting from the introgression time.

### Numerically estimate expected trajectory
Use `./precomp_traj.py` to simulate the expected trajectories. Notice that you can specify several values for the introgression time (which might be significantly faster than making independent simulations for each of those values). Example:
`./precomp_traj.py 0.0006 0.01 10000 1500 1750 2000 2250 > mean_traj.txt`
