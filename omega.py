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

class Omega:
    def __init__(self, pars):
        if not isinstance(pars, list) or len(pars) < 2:
            print("Trajectory should have at least two points.")
            sys.exit(1)
        '''pars.sort(key=lambda interval: interval[0])
        for i in range(len(pars)-1):
            if pars[i][0] == pars[i+1][0]:
                print("All time points should be different.")
                sys.exit(1)'''
        for i in range(len(pars)-1):
            tmp = pars[i][:] + pars[i+1][0:2]
            if pars[i][2] == 0:
                if (pars[i][1] != None and pars[i+1][1] != None and pars[i][1] != pars[i+1][1]) or (pars[i][1] == None and pars[i+1][1] == None):
                    print("Cannot fit trajectory without selection.")
                    sys.exit(1)
                if pars[i+1][0] == None:
                    print("Cannot fit trajectory without selection (time not specified).")
                    sys.exit(1)
                if pars[i][1] == None:
                    pars[i][1] = pars[i+1][1]
                else:
                    pars[i+1][1] = pars[i][1]
                continue
            if tmp.count(None) != 1:
                print("Trajectory cannot be fit, too many or too few free parameters.")
                sys.exit(1)
            if pars[i][0] == None:
                print("T0 cannot be None")
                sys.exit(1)
            elif pars[i][1] == None:
                pars[i][1] = self.time_to_freq(pars[i][0], pars[i+1][0:2]+pars[i][2:])
            elif pars[i][2] == None:
                pars[i][2] = self.fit_select(pars[i:i+2])
            elif pars[i+1][0] == None:
                pars[i+1][0]=self.freq_to_time(pars[i+1][1],pars[i])
            elif pars[i+1][1] == None:
                pars[i+1][1]=self.time_to_freq(pars[i+1][0], pars[i])
        self.intervals = pars
        for i in range(len(self.intervals)-1):
            if self.intervals[i][0] >= self.intervals[i+1][0]:
                raise Exception("trajectory_impossible")
                #print("Time points should be in a sctrictly acsending order.")
                #return(False)
        self.Tp, self.Ta = self.limits()
        self.dT = self.Ta-self.Tp
        self.proportion = self.intervals[-1][1]
                
        
    def limits(self):
        return([self.intervals[0][0], self.intervals[-1][0]])
    
    def time_to_freq(self, t, interval):
        t0 = interval[0]
        omega0 = interval[1]
        s = interval[2]
        p = -s*(t-t0)/2
        k0 = omega0/(1-omega0)
        af = 1.0
        if 1.0-omega0 > 0.0:
            k0 = omega0/(1.0-omega0)
            af = 1.0-1.0/(1.0+k0*exp(p))
        return( af )
    
    def freq_to_time(self, omega, interval):
        if omega == 0.0 or omega == 1.0:
            print("Cannot caclulate time for AF 0.0 or 1.0.")
            sys.exit(1)
        t0 = interval[0]
        omega0 = interval[1]
        if omega0 == 0.0 or omega0 == 1.0:
            print("Cannot caclulate time for AF 0.0 or 1.0.")
            sys.exit(1)
        s = interval[2]
        k0 = omega0/(1-omega0)
        return( t0-2/s*log( 1/k0*omega/(1-omega) ) )
    
    def fit_select(self, pars):
        dT = pars[1][0]-pars[0][0]
        omega0 = pars[0][1]
        omega1 = pars[1][1]
        if omega0 == 0.0 or omega0 == 1.0 or omega1 == 0.0 or omega1 == 1.0:
            print("Cannot caclulate selection coefficient for AF 0.0 or 1.0.")
            sys.exit(1)
        select = -2*(log(omega1/(1 - omega1)) - log(omega0/(1 - omega0)))/dT
        return(select)
    
    def omega(self, t):
        #if t < self.Tp or t > self.Ta:
        #    print("Trajectory is not defined for time", t, "(should be between", self.Tp, "and", self.Ta, ")")
        #    sys.exit(1)
        interval = 0
        while t > self.intervals[interval+1][0] and t < len(self.intervals)-1:
            interval+=1
        return(  self.time_to_freq(t, self.intervals[interval])  )
        
    def Print(self):
        print("time\t\tAF\t\tselection")
        for el in self.intervals[:-1]:
            print("\t\t".join([str(v) for v in el]))
        print("\t\t".join([str(v) for v in self.intervals[-1][:-1]]))