#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

An SIB-model modelling the number of Susceptibles, Infected, and infected with Both malaria strains
individuals over time. 
Inspired by van Baalen (1995) et al. and Alizon (2008) (Multiple infections etc)


Assumptions:

	- No double infections of the same strain occur 
	- S and I1 (the resident strain) are already at equilibrium at the start and static
	- The order of infection does not matter for the double infections

"""
###########################		module import 		###############################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import withinHost_model as inHost

###########################		start parameters	###############################

S0 = 1000 #Susceptibles at time 0
I_1_0 = 50 #Individuals infected with strain 1 (the resident strain)
I_2_0 = 1 #Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both

# number of copies of both strains
n1 = inHost.K1
n2 = inHost.K2
#at equilibrium
n1_12 = (inHost.K1 - inHost.alpha12 * inHost.K2)/(1- inHost.alpha12*inHost.alpha21) 
n2_12 = (inHost.K2 - inHost.alpha21 * inHost.K1)/(1- inHost.alpha12*inHost.alpha21) 

mu = 0.05 #natural death rate

#time
ntimepoints_sys = 20000
time = np.linspace(0,300, ntimepoints_sys)

def eq_sys(y, t, c_delta1, c_delta2, c_beta, h, n1, n2, n1_12, n2_12, mu, S0, I_1):

	#Defining the system of equations 

	beta_1 = inHost.beta_pop(c_beta, n1, h)
	beta_2 = inHost.beta_pop(c_beta, n2, h)
	beta_1_12 = inHost.beta_pop(c_beta, n1_12, h)
	beta_2_12 = inHost.beta_pop(c_beta, n2_12, h)

	delta_1 = inHost.delta_pop(c_delta1, c_delta2, n1, 0)
	delta_2 = inHost.delta_pop(c_delta1, c_delta2, 0, n2)
	delta_12 = inHost.delta_pop(c_delta1, c_delta2, n1_12, n2_12)

	S = S0
	I_1 = I_1_0

	I_2 = y[0]
	I_12 = y[1]

	dI2dt = beta_2*S*I_2 - (mu+delta_2)*I_2 - beta_1*I_1*I_2 

	dI12dt = beta_1*I_1*I_2 + beta_2*I_1*I_2 + beta_2_12*I_1*I_12 - (mu+delta_12)*I_12

	return dI2dt, dI12dt

y0 = (I_2_0, I_12_0)

out = odeint(eq_sys, y0, time, args =(inHost.c_delta1, inHost.c_delta2, inHost.c_beta, 
	inHost.h, n1, n2, n1_12, n2_12, mu, S0, I_1_0))
I_2, I_12 = out.T

############################		Plot system		################################ 

plt.plot(time, I_2, label="Infected with strain 2")
plt.plot(time, I_12, label="Double infection (12)")
plt.legend(loc="best")
plt.xlabel("time")
plt.grid()
plt.show()

plt.semilogy(time, I_2, label="Infected with strain 2")
plt.semilogy(time, I_12, label="Double infection (12)")
#plt.axis([0, 200, 0, 1000])
plt.legend(loc="best")
plt.xlabel("time")
plt.grid()
plt.show()
