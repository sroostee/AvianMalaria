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

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

###########################		start parameters	###############################

S0 = 28 #Susceptibles at time 0
I_1_0 = 1 #Individuals infected with strain 1 (the resident strain)
I_2_0 = 1 #Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both

#transmission rates
beta_1 = 0.02
beta_1_12 = 0.02
beta_12 = 0.02
beta_2 = 0.02
beta_2_12 = 0.02
beta_21 = 0.02
beta_12s = 0.02

#mortality rates caused by parasite
alpha_1 = 0.02
alpha_12 = 0.02
alpha_2 = 0.02

#host natural death rate
mu = 0.02
#host reproduction rate
labda = 0.02
#multiple infection efficiency
epsilon_1 = 0.33
epsilon_2 = 0.33
epsilon_12 = 0.33

#time
ntimepoints = 1000
t = np.linspace(0,200, ntimepoints)

def eq_sys(y, t, beta_1, beta_2, beta_1_12, beta_2_12, beta_12s, 
	alpha_1, alpha_2, alpha_12, mu, labda, epsilon_1, epsilon_2, epsilon_12):

	#Defining the system of equations 

	S = y[0]
	I_1 = y[1]
	I_2 = y[2]
	I_12 = y[3]

	Sdt = (S+I_1+I_2+I_12)*labda - mu*S - beta_1*S*I_1 - epsilon_1*beta_12s*S*I_12 - epsilon_12*beta_12s*S*I_12 - epsilon_2*beta_12s*S*I_2 - beta_2*S*I_2
	I_1dt = beta_1*S*I_1 + epsilon_1*beta_12s*S*I_12 - (mu+alpha_1)*I_1 - beta_12*I_1*I_2 - beta_1_12*I_1*I_12
	I_2dt = beta_2*S*I_2 + epsilon_2*beta_12s*S*I_2 - (mu+alpha_2)*I_2 - beta_21*I_1*I_2 - beta_2_12*I_2*I_12
	I_12dt = beta_12*I_1*I_2 + beta_1_12*I_1*I_12 + beta_21*I_1*I_2 + beta_2_12*I_2*I_12 - (mu+alpha_12)*I_12

	return Sdt, I_1dt, I_2dt, I_12dt

y0 = (S0, I_1_0, I_2_0, I_12_0)

out = odeint(eq_sys, y0, t, args =(beta_1, beta_2, beta_1_12, beta_2_12, beta_12s, 
	alpha_1, alpha_2, alpha_12, mu, labda, epsilon_1, epsilon_2, epsilon_12))
S,I_1,I_2,I_12= out.T

###########################		Plot system		################################ 

plt.plot(t, S, label="Susceptibles")
plt.plot(t, I_1, label="Infected with strain 1")
plt.plot(t, I_2, label="Infected with strain 2")
plt.plot(t, I_12, label="Infected with both")
plt.legend(loc="best")
plt.xlabel("t")
plt.grid()
plt.show()