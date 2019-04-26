#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

An SIB-model modelling the number of Susceptibles, Infected, and infected with Both malaria strains
individuals over time. 
Inspired by van Baalen (1995) et al. and Alizon (2008) (Multiple infections etc)


Assumptions:

	- Within-host dynamics follow Lotka-Volterra competition

"""
###########################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


###########################		start parameters	###############################

n_0 = 1
r0 = 0.1
K = 3

#time
ntimepoints = 1000
t = np.linspace(0,20, ntimepoints)


###########################		growth functions	###############################

def malariaGrowth(n, t, r, K):
	#calculate malaria growth in host
	dndt = r*n*(1-(n/K))
	return dndt

n0 = n_0

out = odeint(malariaGrowth, n0, t, args =(r, K))