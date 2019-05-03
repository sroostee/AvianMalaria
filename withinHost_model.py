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

n1_0 = 1
n2_0 = 1
r1 = 0.2
r2 = 0.5
K1= 10
K2 = 10
alpha12 = 0.1
alpha21 = 0.1

#time
ntimepoints = 1000
t = np.linspace(0,100, ntimepoints)

###########################		growth functions	###############################

def LotkaVolterraCompetition(n, t, r1, r2, K1, K2, alpha12, alpha21):
	#calculate malaria growth in host
	N1 = n[0]
	N2 = n[1]
	dN1dt = r1*N1*(1-((N1+alpha12*N2)/K1))
	dN2dt = r2*N2*(1-((N2+alpha21*N1)/K2))
	return dN1dt, dN2dt

n0 = (n1_0, n2_0)

out = odeint(LotkaVolterraCompetition, n0, t, args =(r1, r2, K1, K2, alpha12, alpha21))

n1, n2 = out.T

###########################		Plot system	over time	################################ 

plt.plot(t, n1, label="Strain 1")
plt.plot(t, n2, label="Strain 2")
plt.legend(loc="best")
plt.xlabel("t")
plt.ylabel("number of copies")
plt.grid()
plt.show()

###########################		Plot system	isoclines	################################ 
plt.plot([0,K2/alpha21],[K2,0], color = "blue", label = "isocline for N1")
plt.plot([0,max(K1, K2/alpha21)+2], [0,0], color = "blue")
plt.plot([0,K1],[K1/alpha12,0], color = "orange", label = "isocline for N2")
plt.plot([0,0], [0,max(K2, K1/alpha12)+2], color = "orange")
plt.axis([-1, max(K1, K2/alpha21)+1, -1, max(K2, K1/alpha12)+1])
plt.legend(loc="best")
plt.grid()
plt.show()

############################  some playing around with direction fields		############


# X = np.arange(-10, 10, 1)
# Y = np.arange(-10, 10, 1)
# U, V = np.meshgrid(X**2, Y**2)

# fig, ax = plt.subplots()
# q = ax.quiver(X, Y, U, V)
# ax.quiverkey(q, X=0.3, Y=1.1, U=10,
#              label='Quiver key, length = 10', labelpos='E')

# plt.show()

# xmax = 4.0
# xmin = -xmax
# D = 20
# ymax = 4.0
# ymin = -ymax
# x = np.linspace(xmin, xmax, D)
# y = np.linspace(ymin, ymax, D)
# X, Y = np.meshgrid(x, y)
# deg = np.arctan(Y**2 - X)
# QP = plt.quiver(X,Y,np.cos(deg),np.sin(deg))
# plt.show()


###########################		Population model 	####################################

c_a = 0.8
c_beta = 0.001
h = 0.5

def a_pop(c_a, n1, n2):
	#calculation on parasite dependent death rate
	# the more copies a host contains the higher the parasite-caused death rate is
	a = c_a*(n1+n2)
	return a 

def beta_pop(c_b, n, h):
	#caculation on parasote dependent transmission rate
	#the more copies a host contains the higher the transmission rate is, but with a max
	beta = (c_b*n)/(1+h*n)
	return beta

###########################		a and beta over time depending on n1 and n2 	#######

a_n1 = a_pop(c_a, n1, 0)
a_n2 = a_pop(c_a, 0, n2)
a_all = a_pop(c_a, n1, n2)
beta_n1 = beta_pop(c_beta, n1, h)
beta_n2 = beta_pop(c_beta, n2, h)

plt.plot(t, a_n1, label = "a n1")
plt.plot(t, a_n2, label = "a n2")
plt.plot(t, a_all, label = "a all")
plt.plot(t, beta_n1, label = "beta_n1")
plt.plot(t, beta_n2, label = "beta_n2")
plt.legend(loc="best")
plt.show()

S0 = 28 #Susceptibles at time 0
I_1_0 = 1 #Individuals infected with strain 1 (the resident strain)
I_2_0 = 1 #Individuals infected with strain 2 (the rare mutant)
I_12_0 = 0 #Individuals infected with both

n1_eq = n1[-1] #number of copies of strain 1 at equilibrium
n2_eq = n2[-1] #number of copies of strain 2 at equilibrium

mu = 0.5 #natural death rate

#time
ntimepoints_sys = 1000
time = np.linspace(0,200, ntimepoints_sys)


def eq_sys(y, t, c_a, c_beta, h, n1, n2, mu, S0, I_1_0):

	#Defining the system of equations 

	beta_1 = beta_pop(c_beta, n1, h)
	beta_2 = beta_pop(c_beta, n2, h)

	a_2 = a_pop(c_a, 0, n2)
	a_12 = a_pop(c_a, n1, n2)

	S = S0
	I_1 = I_1_0

	I_2 = y[0]
	I_12 = y[1]

	dI2dt = beta_2*S*I_2 - (mu+a_2)*I_2 - beta_1*I_1*I_2 

	dI12dt = beta_1*I_1*I_2 + beta_2*I_1*I_2 + beta_2*I_1*I_12 - (mu+a_12)*I_12

	return dI2dt, dI12dt

y0 = (I_2_0, I_12_0)

out = odeint(eq_sys, y0, time, args =(c_a, c_beta, h, n1_eq, n2_eq, mu, S0, I_1_0))
I_2, I_12 = out.T

###########################		Plot system		################################ 

plt.plot(time, I_2, label="Infected with strain 2")
plt.plot(time, I_12, label="Double infection (12)")
plt.axis([0, 200, 0, 10000])
plt.legend(loc="best")
plt.xlabel("t")
plt.grid()
plt.show()