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

User defined functions:
- eq_sys: describing the ODE system of the co-infection model at the start of infection (invasion)

"""
###########################		module import 		###############################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import withinHost_model as inHost

###########################		start parameters	###############################

# number of copies of both strains
n1 = inHost.K1
n2 = inHost.K2
#at equilibrium
n1_12 = (inHost.K1 - inHost.alpha12 * inHost.K2)/(1- inHost.alpha12*inHost.alpha21) 
n2_12 = (inHost.K2 - inHost.alpha21 * inHost.K1)/(1- inHost.alpha12*inHost.alpha21) 

if n1_12 < 0:
	n1_12 = 0

if n2_12 < 0:
	n2_12 = 0

mu = 0.03 #natural death rate
labda = 0.05#birth rate 
l = 1000 #relating to carrying capacity of hosts

delta_1 = inHost.delta_pop(inHost.c_delta1, inHost.c_delta2, n1, 0)
beta_1 = inHost.beta_pop(inHost.c_beta, n1, inHost.h)

S0 = (mu + delta_1)/beta_1 #Susceptibles at time 0
I_1_0 = (-2*labda*(delta_1+mu) - beta_1*l*(delta_1-labda+mu)+
	np.sqrt(beta_1)*np.sqrt(l) * np.sqrt(4*delta_1*labda* (delta_1+mu) + beta_1*l*(delta_1-labda+mu)**2))/(2*beta_1*labda)
#Individuals infected with strain 1 (the resident strain)
I_2_0 = 1#Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both

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

	dI2dt = beta_2*S*I_2 + beta_2_12*S*I_12 - (mu+delta_2)*I_2 - beta_1*I_1*I_2 

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
