#!/usr/bin/env python3

#    Copyright (c) 2019 Russel Corbett-Detig, Vladimir Shchur (vlshchur@gmail.com)
#
#    This file is part of DAIM.
#
#    DAIM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DAIM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DAIM.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import sys
import argparse


parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('proportion', type=float,
                    help='admixture proportion')
parser.add_argument('selection', type=float,
                    help='selection coefficient (for homozygous diploid individual with fitness 1+s, which translates into the fitness 1+s/2 for haploid populations)')
parser.add_argument('Ne', type=int,
                    help='haploid effective population size')
parser.add_argument('generations', type=int, nargs='+',
                    help='number of generations (can have multiple values)')

parser.add_argument('-at', action='store_true',
                    help='Output all trajectories, if not specified, only expected trajectories will be output')
parser.add_argument('--reps', '-r', type=int, default=100000,
                    help='number of repetitions')
                    

                    
clargs = parser.parse_args()
if isinstance(clargs.at, list):
    clargs.at = clargs.at[0]
if isinstance(clargs.reps, list):
    clargs.reps = clargs.reps[0]

clargs.selection /= 2.0

def SimulateTrajectory(f, s, n, g, output_tr=False):
    na = n * f
    nA = n * ( 1 - f )

    result = [f]
    fix_freq = None

    for gen in range(0,g+1) :

	### total fitness for each genotype
        na_fitness = na * ( 1 + s )
        nA_fitness = nA

        draw_p = float( na_fitness ) / float( na_fitness + nA_fitness )
        if draw_p == 0.0 or draw_p == 1.0:
            fix_freq = draw_p
            break

        new_na = np.random.binomial(n,draw_p)

        result.append(float(new_na)/float(n))

        na = new_na
        nA = n - na
    while len(result) < g+1:
        result.append(fix_freq)
    if output_tr:
        print( "TR\t" + "\t".join(result) )
    return result
    
    
def CalculateMeanTrajectories(trajs, gens):
    mean_traj = { g: [0.0 for _ in range(g+1)] for g in gens }
    counters = {g: 0 for g in gens}
    for tr in trajs:
        for g in gens:
            if tr[g] != 0.0:
                counters[g] += 1
                for i in range(g+1):
                    mean_traj[g][i] += tr[i]
    sys.stderr.write("Number of sampled trajectories:\n")
    for g in gens:
        for i in range( len(mean_traj[g]) ):
            mean_traj[g][i] /= counters[g]
        sys.stderr.write(str(g) + " generations: " + str(counters[g]) + " trajectories\n")
    return mean_traj

## results vectors
results = []
max_gen = max(clargs.generations)

for sim in range(clargs.reps) :
    result = SimulateTrajectory(clargs.proportion, clargs.selection, clargs.Ne, max_gen, clargs.at)
    results.append(result)
mean_traj = CalculateMeanTrajectories(results, clargs.generations)


for g, tr in mean_traj.items():
    print( "MT\t" + "\t".join( [str(v) for v in tr] ) )

