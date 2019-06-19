#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

Modelling the dynamics of within-host competition of malaria parasites in avian malaria.
The dynamics follow Lotka-Volterra competition of two competitors which has the form of:
dNidt = ri*Ni*(1-((Ni+alphaij*Nj)/Ki))
dNjdt = rj*Nj*(1-((Nj+alphaji*Ni)/Kj))

Assumptions:

	- Within-host dynamics follow Lotka-Volterra competition
	- Only two competitors are present

User defined functions:
- LotkaVolterraCompetition: function for the ODE of Lotka-Volterra competition dynamics
- delta_pop: calculation of delta for the host population through values of n1 and n2
- beta_pop: calculation of beta for the host population through the values of n1 and n2

"""
###########################		module import 		###############################

import numpy as np
from scipy.integrate import odeint
import matplotlib
import matplotlib.pyplot as plt

###########################		start parameters	###############################

n1_0 = 1
n2_0 = 1
r1 = .5
r2 = .5
K1= 5
K2 = 5
alpha12 = .5
alpha21 = .5

n0 = (n1_0, n2_0)

#time
ntimepoints = 100
t = np.linspace(0,100, ntimepoints)

############ parameters for delta and beta in host model
c_delta1 = 0.005
c_delta2 = 0.005
c_beta = 0.0001
h = 0.05

###########################		growth functions	###############################

def LotkaVolterraCompetition(n, t, r1, r2, K1, K2, alpha12, alpha21):
	#calculate malaria growth in host
	N1 = n[0]
	N2 = n[1]
	dN1dt = r1*N1*(1-((N1+alpha12*N2)/K1))
	dN2dt = r2*N2*(1-((N2+alpha21*N1)/K2))
	return dN1dt, dN2dt

###########	Effect of within-host dynamics on delta and beta in the population model 	###########

def delta_pop(c_delta1, c_delta2, n1, n2):
	#calculation on parasite dependent death rate
	# the more copies a host contains the higher the parasite-caused death rate is
	delta = c_delta1*n1+c_delta2*n2
	return delta 

def beta_pop(c_b, n, h):
	#caculation on parasote dependent transmission rate
	#the more copies a host contains the higher the transmission rate is, but with a max
	beta = (c_b*n)/(1+(h*n))
	return beta


# only do this if it is the main plot
if __name__ == "__main__":

	######################		run ODE system		##################################################

	out = odeint(LotkaVolterraCompetition, n0, t, args =(r1, r2, K1, K2, alpha12, alpha21))

	n1, n2 = out.T

	n0_single1 = (1,0)
	n0_single2 = (0,1)

	out_singlen1 = odeint(LotkaVolterraCompetition, n0_single1, t, args =(r1, r2, K1, K2, alpha12, alpha21))
	out_singlen2 = odeint(LotkaVolterraCompetition, n0_single2, t, args =(r1, r2, K1, K2, alpha12, alpha21))

	n1_single, n2_non = out_singlen1.T
	n1_non, n2_single = out_singlen2.T

	############################	increase font size for plotting		#######################

	font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

	matplotlib.rc('font', **font)  # pass in the font dict as kwargs


	###########################		delta and beta over time depending on n1 and n2 	#######

	delta_n1 = delta_pop(c_delta1, c_delta2, n1_single, 0)
	delta_n2 = delta_pop(c_delta1, c_delta2, 0, n2_single)
	delta_co = delta_pop(c_delta1, c_delta2, n1, n2)
	beta_n1 = beta_pop(c_beta, n1_single, h)
	beta_n2 = beta_pop(c_beta, n2_single, h)
	beta_n1_co = beta_pop(c_beta, n1, h)
	beta_n2_co = beta_pop(c_beta, n2, h)

###########################		Plot within-host system	over time	################################ 

	plt.plot(t, n1_single, 's', color = 'blue' , label="strain 1")
	plt.plot(t, n2_single, ':', color = 'orange' , label="strain 2")
	plt.plot(t, n1, 'P', color = "blue", label="strain 1 co")
	plt.plot(t, n2, '-', color = 'orange' , label="strain 2 co")
	plt.legend(loc="best")
	plt.xlabel("t")
	plt.ylabel("number of copies")
	plt.grid()
	plt.show()

	###########################		Plot system	isoclines	########################################## 
	plt.plot([0,K2/alpha21],[K2,0], color = "blue", label = "n1")
	plt.plot([0,max(K1, K2/alpha21)+2], [0,0], color = "blue")
	plt.plot([0,K1],[K1/alpha12,0], color = "orange", label = "n2")
	plt.plot([0,0], [0,max(K2, K1/alpha12)+2], color = "orange")
	plt.axis([-1, max(K1, K2/alpha21)+1, -1, max(K2, K1/alpha12)+1])
	plt.legend(loc="best")
	plt.xlabel("n1")
	plt.ylabel("n2")
	plt.grid()
	plt.show()

	##########################		Plot delta values over time ######################################

	plt.plot(t, delta_n1, 's', color = "blue", label = r'$\delta$ n1')
	plt.plot(t, delta_n2, ':', color = "orange", label = r'$\delta$ n2')
	plt.plot(t, delta_co, '-.', color = "deeppink", label = r'$\delta$ co')
	plt.legend(loc="best")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel("t")
	plt.ylabel(r'$\delta$')
	plt.grid()
	plt.show()

	###########################3	Plot beta values over time ########################################

	plt.plot(t, beta_n1, 's',  color = "blue", label = r'$\beta$ n1')
	plt.plot(t, beta_n2,':',  color = "orange", label = r'$\beta$ n2')	
	plt.plot(t, beta_n1_co, 'P',  color = "blue", label = r'$\beta$ n1 co')
	plt.plot(t, beta_n2_co, '-', color = "orange", label = r'$\beta$ n2 co')
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	plt.xlabel("t")
	plt.legend(loc="best")
	plt.ylabel(r'$\beta$')
	plt.grid()
	plt.show()