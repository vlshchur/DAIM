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
from omega import Omega
from ode import ExpectedTractLength



Ne = 10000#Effective population size
trajectory_parameters = [[0, None, 0.01], [500, 0.1, 0]]#Allele trajectory



omega = Omega(trajectory_parameters)
omega.Print()
if False:
    plt.plot(range(omega.Tp,omega.Ta), [omega.omega(t) for t in range(omega.Tp,omega.Ta)])
    plt.show()

tr_len = ExpectedTractLength(omega, Ne)
print("tract length is ", tr_len)