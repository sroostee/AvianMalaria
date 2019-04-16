#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

A simple SID-model modelling the number of Susceptibles, Infected and Double infected
individuals over time. 
Following van Baalen (1995) et al. and Alizon (2008) (Multiple infections etc)


Assumption:

No double infection by the mutant strain as this is a rare strain
"""
###########################


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#start parameters

S0 = 24
I_1_0 = 5
I_2_0 = 1
D_11_0 = 0
D_12_0 = 0
D_21_0 = 0

#transmission rates
beta_1 = 0.02
beta_1_11 = 0.02
beta_1_12 = 0.02
beta_1_21 = 0.02
beta_2 = 0.02
beta_2_12 = 0.02
beta_2_21 = 0.02

#mortality rates caused by parasite
alpha_1 = 0.01
alpha_11 = 0.01
alpha_12 = 0.02
alpha_2 = 0.015
alpha_21 = 0.03

#host natural death rate
mu = 0.02
#host reproduction rate
rho = 0.05
#multiple infection efficiency
epsilon = 0.2

#time
ntimepoints = 1000
t = np.linspace(0,20, ntimepoints)

def ode_multiple(y, t, beta_1, beta_2, beta_1_11, beta_1_12, beta_1_21, beta_2_12, 
	beta_2_21, alpha_1, alpha_2, alpha_11, alpha_12, alpha_21, mu, rho, epsilon):
	S = y[0]
	I_1 = y[1]
	I_2 = y[2]
	D_11 = y[3]
	D_12 = y[4]
	D_21 =y[5]
	#force of infection rates for strain 1 and 2 
	labda_1 = beta_1*I_1 + beta_1_11*D_11 + beta_1_12*D_12 + epsilon*beta_1_21*D_21
	labda_2 = beta_2*I_2 + epsilon*beta_2_12*D_12 + beta_2_21*D_21

	Sdt = (S+I_1+I_2+D_11+D_12+D_21)*rho - (labda_1+labda_2+mu)*S
	I_1dt = labda_1*S - (alpha_1+labda_1+labda_2+mu)*I_1
	I_2dt = labda_2*S - (alpha_2+labda_1+labda_2+mu)*I_2
	D_11dt = labda_1*I_1 - (alpha_11+mu)*D_11
	D_12dt = labda_2*I_1 - (alpha_12+mu)*D_12
	D_21dt = labda_1*I_2 - (alpha_21+mu)*D_21

	return Sdt, I_1dt, I_2dt, D_11dt, D_12dt, D_21dt

y0 = (S0, I_1_0, I_2_0, D_11_0, D_12_0, D_21_0)


out = odeint(ode_multiple, y0, t, args =(beta_1, beta_2, beta_1_11, beta_1_12, beta_1_21, beta_2_12, 
	beta_2_21, alpha_1, alpha_2, alpha_11, alpha_12, alpha_21, mu, rho, epsilon))
S,I_1,I_2,D_11,D_12,D_21= out.T

plt.plot(t, S, label="Susceptibles")
plt.plot(t, I_1, label="Infected with strain 1")
plt.plot(t, I_2, label="Infected with strain 2")
plt.plot(t, D_11, label="Double infection (11)")
plt.plot(t, D_12, label="Double infection (12)")
plt.plot(t, D_21, label="Double infection (21)")
plt.legend(loc="best")
plt.xlabel("t")
plt.grid()
plt.show()