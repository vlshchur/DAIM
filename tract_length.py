#!/usr/bin/env python3

#    Copyright (c) 2019 Vladimir Shchur (vlshchur@gmail.com)
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

import sys
import os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import matplotlib.pyplot as plt
from omega import Omega_logit,Omega_precomp
from ode import ExpectedTractLength

def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('mode', choices=['logit', 'precomp'],
                    help='trajectory mode')
parser.add_argument('input_file', nargs='*',
                    help='input file with the trajectory for precomputed mode')

parser.add_argument('--point', '-p', nargs=3, type=str, action = 'append',
                    help='points for logistic function mode: [time] [AF] [selectio]')#-mi [npop:1/2] [migStart] [migEnd] [init val] [var:0/1]

parser.add_argument('-Ne', type=int, default=10000,
                    help='haploid effective population size')

parser.add_argument('--pdf', action='store_true',
                    help='Output probability density function')

#parser.add_argument('-at', action='store_true',
#                    help='Output all trajectories, if not specified, only expected trajectories will be output')                    

                    
cl = parser.parse_args()
if isinstance(cl.Ne, list):
    cl.Ne = cl.Ne[0]

trajectory_parameters = None
omega = None
lens = []
if cl.mode == 'logit':
    if cl.point is None or len(cl.point) < 2:
        print("Logistic mode needs at least two points to set the trajectory.")
        sys.exit(1)
    trajectory_parameters = [[None if v=="None" else float(v) for v in p] for p in cl.point]
    try:
        omega = Omega_logit(trajectory_parameters)
        tr_len = ExpectedTractLength(omega, cl.Ne)
        lens.append([omega.dT, tr_len])
    except ValueError:
        print("Cannot initialise AF trajectory, allele goes to fixation.")
        sys.exit(1)
elif cl.mode == 'precomp':
    if len(cl.input_file) == 0:
        print("Please specify at least one input file with trajectories.")
        sys.exit(1)
    for fn in cl.input_file:
        with open(fn) as infile:
            for line in infile:
                line = line.split("\t")
                if line[0] == "MT":
                    omega = Omega_precomp([float(v) for v in line[1:]])
                    tr_len = ExpectedTractLength(omega, cl.Ne)
                    lens.append([omega.dT, tr_len])
else:
    print("Unknown mode...")
    sys.exit(1)


print("time\tlength\tsd")
for v in lens:
    print(v[0], "\t", v[1][0], "\t", v[1][1])
    if cl.pdf:
        print_err("#HEADER", v[0], sep="\t")
        for prob, point in zip(v[1][2], v[1][3]):
            print_err(point, prob, sep="\t")

sys.exit(0)


















Ne = 10000#Effective population size
times = [1500,1750,2000,2250]
times = [1500]
lens = []
trajectory_parameters = [[0, None, 0.01], [1500, 0.0006, 0]]#Allele trajectory

fn_input = "data/denisova_means.txt"

for maxT in times:
    try:
#        omega = Omega1(0.01,0.0006,maxT,10000)
#        omega = Omega(trajectory_parameters)
        omega = Omega_matrix(fn_input, maxT)
    except ValueError:
        print("Cannot initialise AF trajectory, allele goes to fixation.")
        print("time\tlength")
        for t, l in zip(times, len):
            print(t, "\t", l)
        sys.exit(1)
    omega.Print()
    if False:
    #    plt.plot(range(omega.Tp,omega.Ta), [omega.omega(t) for t in range(omega.Tp,omega.Ta)])
        plt.plot(range(0,maxT), [omega.omega(t) for t in range(0,maxT)])
#        plt.plot(range(0,maxT), [omega.omega1(t) for t in range(0,maxT)])
        plt.plot(range(0,maxT), [omega.omega2(t) for t in range(0,maxT)])
#        plt.plot(range(0,maxT), [omega.omega3(t) for t in range(0,maxT)])
        plt.show()
        sys.exit(0)

    #print("Calculating...")

    tr_len = ExpectedTractLength(omega, Ne)
    lens.append(tr_len)
    #print("tract length is ", tr_len)
print("time\tlength")
for t, l in zip(times, lens):
    print(t, "\t", l)

sys.exit(0)

try:
    omega = Omega(trajectory_parameters)
except ValueError:
    print("Cannot initialise AF trajectory, allele goes to fixation.")
    sys.exit(1)
    
print("Calculating...")

tr_len = ExpectedTractLength(omega, Ne)

print("tract length is ", tr_len)