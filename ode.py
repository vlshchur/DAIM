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
#    along with MiSTI.  If not, see <https://www.gnu.org/licenses/>.

import sys
import os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import numpy as np
import math
from math import (exp,log)
from scipy import (linalg as la,integrate)
import time
from numpy import (dot,identity)
import matplotlib.pyplot as plt
#from omega import Omega

class precision:
    def __init__(self):
        self.r2_c = 1.0
        self.range_c = 5.0
        self.discr = 1000
        

class ODE_parameters:
    def __init__(self, Ne, r1, r2=0.000001):
        self.Ne = Ne
        self.r1 = r1
        self.r2 = r2

def ODE(y, t, omega, ode_pars):
    lam = 1/ode_pars.Ne
    r1 = ode_pars.r1
    r2 = ode_pars.r2
    om = omega(t)
    M = [
            [-r1*(1 - om) - r2, lam, 0, 0, 0, r1*om],
            [r2*om, -lam - (2*r1 + r2)*(1 - om), (r1 + r2)*om, r1*om, 0, 0 ],
            [r2*(1 - om), (r1 + r2)*(1 - om), -r1 - r2*om, 0, r1*om, 0],
            [0, r1*(1 - om), 0, -r1 - r2*(1 - om), (r1 + r2)*om, r2*om],
            [0, 0, r1*(1 - om), (r1 + r2)*(1 - om), -lam - (2*r1 + r2)*om, r2*(1 - om)],
            [r1*(1 - om), 0, 0, 0, lam, -r1*om - r2]
        ]
    return( dot(M, y) )

def ExpectedTractLength(omega, Ne):
    prec = precision()
    rRange = prec.range_c*NeutralExpectation( Ne, omega.dT, omega.proportion )
    '''print("rRange = ", rRange)
    print("Integration interval Tp = ", Tp, "\tTa = ", Ta)
    print("Frequencies omega_p = ", time_to_freq(Tp, s), "\tomega_a = ", time_to_freq(Ta, s))'''
    dr = rRange/prec.discr
    ode_pars = ODE_parameters(Ne, None, dr*prec.r2_c)#selection coefficient, coalescent rate (1/effective population size), r1, r2, where r1 is the distance from the selected site to the first loci and r2 is the distance between two loci
    transition_rates = []
    p0 = [1, 0, 0, 0, 0, 0]
    for d in range(prec.discr):
        ode_pars.r1 = dr*d
        sol = integrate.odeint(lambda y, t: ODE(y, t, omega.omega, ode_pars), p0, omega.limits())
        transition_rates.append( sol[1][2]/(sol[1][0]+sol[1][1]+sol[1][2])/ode_pars.r2 )
    cumul_rate = 0
    pdf = []
    for d in range(prec.discr):
       pdf.append(exp(-cumul_rate*dr)*transition_rates[d])
       cumul_rate += transition_rates[d]
    norm = sum(pdf)
    expected_tr_len = 0
    for d in range(prec.discr):
        expected_tr_len += pdf[d]*(d*dr)
    return( 2*expected_tr_len/norm )

def time_to_freq(t, s):
    return(  1 - 1/(1 + exp(-s*t/2))  )

def NeutralExpectation(Ne, T, omega1):
    return(  2/(2*Ne*(1 - omega1)*(1 - exp(-T/2/Ne)))  )

def freq_to_time(omega, s):
    return(  -2/s*log(omega/(1 - omega))  )