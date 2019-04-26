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
###########################

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

###########################		start parameters	###############################

S0 = 28 #Susceptibles at time 0
I_1_0 = 1 #Individuals infected with strain 1 (the resident strain)
I_2_0 = 1 #Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both

#transmission rates
#beta_1 = 0.02
beta_12 = 0.02
beta_1_12 = 0.02
#beta_1_21 = 0.02
beta_2 = 0.02
beta_2_12 = 0.02
#beta_2_21 = 0.02
beta_12s = 0.02
beta_21 = 0.02

#mortality rates caused by parasite
#alpha_1 = 0.02
#alpha_11 = 0.02
alpha_12 = 0.02
alpha_2 = 0.02
#alpha_21 = 0.02

#host natural death rate
mu = 0.02
#host reproduction rate
#labda = 0.02
#multiple infection efficiency: epsilon_1 + epsilon_2 + epsilon_12 = 1
#epsilon_1 = 0.33
epsilon_2 = 0.33
epsilon_12 = 0.33

#time
ntimepoints = 1000
t = np.linspace(0,20, ntimepoints)

def eq_sys(y, t, beta_2, beta_12s, beta_21, beta_2_12, epsilon_2, mu, alpha_2, 
	epsilon_12, beta_12, beta_1_12, alpha_12, S0, I_1_0):

	#Defining the system of equations 

	S = S0
	I_1 = I_1_0
	#S = y[0]
	#I_1 = y[1]
	I_2 = y[0]
	I_12 = y[1]

	#dSdt = S0
	#dI1dt = I_1_0

	dI2dt = beta_2*S*I_2 + epsilon_2*beta_12s*S*I_12 
	- (mu+alpha_2)*I_2 - beta_21*I_1*I_2 - beta_2_12*I_2*I_12

	dI12dt = epsilon_12*beta_12s*S*I_12 + beta_21*I_1*I_2 + beta_2_12*I_2*I_12 
	+ beta_12*I_1*I_2 + beta_1_12*I_1*I_12 - (mu+alpha_12)*I_12

	return dI2dt, dI12dt

y0 = (I_2_0, I_12_0)

# out = odeint(eq_sys, y0, t, args =(beta_2, beta_12s, beta_21, beta_2_12, epsilon_2, mu, alpha_2, 
# 	epsilon_12, beta_12, beta_1_12, alpha_12, S0, I_1_0))
# I_2, I_12 = out.T


###########################		Plot system		################################ 
#plt.plot(t, S, label="Susceptibles")
#plt.plot(t, I_1, label="Infected with strain 1")
# plt.plot(t, I_2, label="Infected with strain 2")
# plt.plot(t, I_12, label="Double infection (12)")
# plt.axis([0, 15, 0, 500])
# plt.legend(loc="best")
# plt.xlabel("t")
# plt.grid()
# plt.show()

###########################		Calculate R0	################################

NGM = np.matrix([[]])

transmissions = np.matrix([[beta_2*S0,epsilon_2*beta_12s*S0],[beta_21*I_1_0,beta_2_12*I_2_0]])

transitions = np.matrix([[mu+alpha_2, 0], [0, mu+alpha_12]])  

K = -transmissions * transitions.I

print(transmissions)
print(transitions)
print(K)