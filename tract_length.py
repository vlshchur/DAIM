#!/usr/bin/env python3

#    Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)
#
#    This file is part of AITL.
#
#    AITL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AITL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AITL.  If not, see <https://www.gnu.org/licenses/>.

import sys
import os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import matplotlib.pyplot as plt
from omega import Omega,Omega1,Omega_matrix
from ode import ExpectedTractLength



Ne = 10000#Effective population size
times = [1500,1750,2000,2250]
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
        plt.show()

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