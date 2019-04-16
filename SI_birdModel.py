#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

A simple SIR-model modelling the number of Susceptibles, Infected  individuals
over time. 


"""
###########################


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N_0 = 1000
I1_0 = 1
I2_0 = 1
I12_0 = 0
S0 = N_0 - I1_0 - I2_0 - I12_0

#contact rates
beta_1 = 0.08 #contact+infection rate of strain 1
beta_2 = 0.2 #contact+infection rate of strain 2
beta_12 = 0.1 #contact+infection rate of already infected strain 1 to strain 2
beta_21 = 0.1 #contact+infection rate of already infected strain 2 to strain 1

#virulence rates
alpha_1 = 0.9
alpha_2 = 0.01
alpha_12 = 0.01

#vital dynamics
labda = 0.07 #birth rate
mu = 0.07 #death rate

#time
ntimepoints = 1000
t = np.linspace(0,300, ntimepoints)

# def deriv(y, t, beta_1, beta_2, labda, mu, alpha_1, alpha_2):
# 	S = y[0]
# 	I_1 = y[1]
# 	I_2 = y[2]
# 	I_12 = y[3]
#	N = S+I_1+I_2+I_12
# 	dSdt = labda*N -(beta_1*S*I_1)/N -(beta_2*S*I_2)/N - mu*S
# 	dI_1dt = (beta_1*S*I_1)/N - (mu+alpha_1)*I_1
# 	dI_2dt = (beta_2*S*I_2)/N - (mu+alpha_2)*I_2
# 	dI12_dt = (beta_12*I_1*I_2)/N + (beta_21*I_2*I_1)/N - (mu+alpha_12)*I_12
# 	#return dSdt, dI_1dt, dI_2dt, dI12_dt

def deriv(y, t, beta_1, labda, mu, alpha_1):
	S = y[0]
	I_1 = y[1]
	N = S + I_1
	R_0 = beta_1*
	dSdt = labda*N -(beta_1*S*I_1)/N - mu*S
	dI_1dt = (beta_1*S*I_1)/N - (mu+alpha_1)*I_1
	R_0dt = (beta_1*(I_1/N))/(mu+alpha_1)
	return dSdt, dI_1dt, R_0dt
	#return dSdt, dI_1dt, dI_2dt, dI12_dt

y0 = (S0, I1_0, 0)
#y0 = (S0, I1_0, I2_0, I12_0)
#out = odeint(deriv, y0, t, args =(beta_1, beta_2, labda, mu, alpha_1, alpha_2))
#S,I_1,I_2,I_12= out.T
out = odeint(deriv, y0, t, args =(beta_1,labda, mu, alpha_1))
S,I_1,R_0= out.T

plt.plot(t, S, label="Susceptibles")
plt.plot(t, I_1, label="Infected with strain 1")
#plt.plot(t, I_2, label="Infected with strain 2")
#plt.plot(t, I_12, label="Infected with both")
plt.plot(t, R_0, label="Reproduction ratio")
plt.legend(loc="best")
plt.xlabel("t")
plt.grid()
plt.show()

############################################

#R_0 = (beta_1*S0)/(mu+alpha_1)


#np.log(s_inf - 1) - np.log(s_inf)