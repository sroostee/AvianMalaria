#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

An SIB-model modelling the number of Susceptibles, Infected, and infected with Both malaria strains
individuals over time. 
Inspired by van Baalen (1995) et al. and Alizon (2008) (Multiple infections etc)

"""
###########################


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import withinHost_model as inHost

###########################		start parameters	###############################

S_0 = 299 #Susceptibles at time 0
I_1_0 = 1 #Individuals infected with strain 1 (the resident strain)
I_2_0 = 1 #Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both

n1_eq = inHost.n1[-1] #number of copies of strain 1 at equilibrium
n2_eq = inHost.n2[-1] #number of copies of strain 2 at equilibrium

mu = 0.05 #natural death rate
labda = 0.9#birth rate

#time
ntimepoints_sys = 1000
time = np.linspace(0,40, ntimepoints_sys)

def eq_sys(y, t, c_delta, c_beta, h, n1, n2, mu, labda):

	#Defining the system of equations 

	beta_1 = inHost.beta_pop(c_beta, n1, h)
	beta_2 = inHost.beta_pop(c_beta, n2, h)

	delta_1 = inHost.delta_pop(c_delta, n1, 0)
	delta_2 = inHost.delta_pop(c_delta, 0, n2)
	delta_12 = inHost.delta_pop(c_delta, n1, n2)

	S = y[0]
	I_1 = y[1]
	I_2 = y[2]
	I_12 = y[3]

	dSdt = labda * (S+I_1+I_2+I_12) - mu*S -beta_1*S*I_1 - beta_2*S*I_2

	dI1dt = beta_1*S*I_1 - (mu+delta_1)*I_1 - beta_2*I_1*I_2 - beta_2*I_1*I_12

	dI2dt = beta_2*S*I_2 - (mu+delta_2)*I_2 - beta_1*I_1*I_2 

	dI12dt = beta_1*I_1*I_2 + beta_2*I_1*I_2 + beta_2*I_1*I_12 - (mu+delta_12)*I_12

	return dSdt, dI1dt, dI2dt, dI12dt

y0 = (S_0, I_1_0, I_2_0, I_12_0)

out = odeint(eq_sys, y0, time, args =(inHost.c_delta, inHost.c_beta, inHost.h, n1_eq, n2_eq, mu, labda))
S, I_1, I_2, I_12 = out.T

###########################		Plot system		################################ 
#only plot if this is the main script
if __name__ == "__main__":

	plt.plot(time, S, label="Susceptibles")
	plt.plot(time, I_1, label="Infected with strain 1")
	plt.plot(time, I_2, label="Infected with strain 2")
	plt.plot(time, I_12, label="Infected with both")
	plt.legend(loc="best")
	plt.xlabel("t")
	plt.grid()
	plt.show()

