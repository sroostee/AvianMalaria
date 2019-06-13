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

# number of copies of both strains
n1 = inHost.K1
n2 = inHost.K2
#at equilibrium
n1_12 = (inHost.K1 - inHost.alpha12 * inHost.K2)/(1- inHost.alpha12*inHost.alpha21) 
n2_12 = (inHost.K2 - inHost.alpha21 * inHost.K1)/(1- inHost.alpha12*inHost.alpha21) 

if n1_12 < 0:
	n1_12 = 0
	n2_12 = inHost.K2

if n2_12 < 0:
	n2_12 = 0
	n1_12 = inHost.K1

mu = 0.03 #natural death rate
labda = 0.05#birth rate 
l = 1000 #relating to carrying capacity of host population

delta_1 = inHost.delta_pop(inHost.c_delta1, inHost.c_delta2, n1, 0)
beta_1 = inHost.beta_pop(inHost.c_beta, n1, inHost.h)

S_0 = (mu + delta_1)/beta_1 #Susceptibles at time 0
I_1_0 = (-2*labda*(delta_1+mu) - beta_1*l*(delta_1-labda+mu)+
	np.sqrt(beta_1)*np.sqrt(l) * np.sqrt(4*delta_1*labda* (delta_1+mu) + beta_1*l*(delta_1-labda+mu)**2))/(2*beta_1*labda)
#Individuals infected with strain 1 (the resident strain)
I_2_0 = 1#Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both
DI_0 = 0 #deceased hosts through infection
D_0 = 0 #deceased hosts total

#time
ntimepoints_sys = 20000
time = np.linspace(0,2000, ntimepoints_sys)

def eq_sys(y, t, c_delta1, c_delta2, c_beta, h, n1, n2, n1_12, n2_12, mu, labda, l):

	#Defining the system of equations 

	beta_1 = inHost.beta_pop(c_beta, n1, h)
	beta_2 = inHost.beta_pop(c_beta, n2, h)
	beta_1_12 = inHost.beta_pop(c_beta, n1_12, h)
	beta_2_12 = inHost.beta_pop(c_beta, n2_12, h)

	delta_1 = inHost.delta_pop(c_delta1, c_delta2, n1, 0)
	delta_2 = inHost.delta_pop(c_delta1, c_delta2, 0, n2)
	delta_12 = inHost.delta_pop(c_delta1, c_delta2, n1_12, n2_12)

	S = y[0]
	I_1 = y[1]
	I_2 = y[2]
	I_12 = y[3]
	DI = y[4]
	D = y[5]

	dSdt = labda*(S+I_1+I_2+I_12) *(1-(S+I_1+I_2+I_12)/l)  - mu*S -beta_1*S*I_1 - beta_2*S*I_2 - beta_1_12*S*I_12 - beta_2_12*S*I_12

	dI1dt = beta_1*S*I_1 + beta_1_12*S*I_12 - (mu+delta_1)*I_1 - beta_2*I_1*I_2 - beta_2_12*I_1*I_12

	dI2dt = beta_2*S*I_2 + beta_2_12*S*I_12 - (mu+delta_2)*I_2 - beta_1*I_1*I_2 - beta_1_12*I_2*I_12

	dI12dt = beta_1*I_1*I_2 + beta_2*I_1*I_2 + beta_2_12*I_1*I_12 + beta_1_12*I_2*I_12 - (mu+delta_12)*I_12

	dDIdt = delta_1*I_1 + delta_2*I_2 + delta_12*I_12

	dDdt = mu*S + (mu+delta_1)*I_1 + (mu+delta_2)*I_2 + (mu+delta_12)*I_12

	return dSdt, dI1dt, dI2dt, dI12dt, dDIdt, dDdt

y0 = (S_0, I_1_0, I_2_0, I_12_0, DI_0, D_0)

out = odeint(eq_sys, y0, time, args =(inHost.c_delta1, inHost.c_delta2, inHost.c_beta, 
	inHost.h, n1, n2, n1_12, n2_12, mu, labda, l))
S, I_1, I_2, I_12, DI, D = out.T

###########################		Plot system		################################ 
#only plot if this is the main script
if __name__ == "__main__":

	plt.plot(time, S, label="Susceptibles")
	plt.plot(time, I_1, label="Infected with strain 1")
	plt.plot(time, I_2, label="Infected with strain 2")
	plt.plot(time, I_12, label="Infected with both")
	plt.axis([0,2000,0,300])
	plt.legend(loc="best")
	plt.xlabel("t")
	plt.grid()
	plt.show()

	plt.plot(time, D, label = "Deceased hosts")
	plt.plot(time, DI, label = "Deceased hosts though infection")
	plt.legend(loc="best")
	plt.grid()
	plt.show()

