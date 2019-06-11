#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

Script written to visualise the effect of the Lotka-Volterra competition parameters on the 
number of double infected hosts. 

"""
###########################

import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import seaborn as sns
import pandas as pd
import withinHost_model as inHost
import host_model as host

######################### generate alpha values ############################
alpha_ij = np.arange(0.1,1.0,0.1)
alpha_ji = np.arange(0.1,1.0,0.1)

i_1 = []
i_2 = []
i_ij = []
dead = []
a_ij = []
a_ji = []

for a_i in alpha_ij:
	for a_j in alpha_ji:
		a_ij.append(a_i)
		a_ji.append(a_j)

		#calculate copies of each strain in the host
		if a_i == 1 and a_j == 1:
			#special case where intra- and interspecies competition is the same
			n_out = odeint(inHost.LotkaVolterraCompetition, inHost.n0, inHost.t, args =(inHost.r1, inHost.r2, inHost.K1, inHost.K2, a_i, a_j))
			n1_12, n2_12 = n_out.T
			n1_12 = n1_12[-1]
			n2_12 = n2_12[-1]
		else:
			#at equilibrium
			n1_12 = (inHost.K1 - a_i * inHost.K2)/(1- a_i*a_j) 
			n2_12 = (inHost.K2 - a_j * inHost.K1)/(1- a_i*a_j)

		n1 = inHost.K1
		n2 = inHost.K2 
		#calculate the number of hosts with a double infection
		host_out = odeint(host.eq_sys, host.y0, host.time, args =(inHost.c_delta1, inHost.c_delta2, inHost.c_beta, inHost.h, n1, n2, n1_12, n2_12, host.mu, host.labda, host.l))
		S, I_1, I_2, I_12, dead_host = host_out.T
		I_1 = I_1[-1]
		I_2 = I_2[-1]
		I_12 = I_12[-1]
		dead_host = dead_host[-1]
		i_1.append(I_1)
		i_2.append(I_2)
		i_ij.append(I_12)
		dead.append(dead_host)


I1_alpha_df = pd.DataFrame(dict(alphaij = a_ij, alphaji = a_ji, I1 = i_1))
I1_to_alphas = I1_alpha_df.pivot("alphaij", "alphaji", "I1")
# print(I1_alpha_df)
# print(I1_to_alphas)

I2_alpha_df = pd.DataFrame(dict(alphaij = a_ij, alphaji = a_ji, I2 = i_2))
I2_to_alphas = I2_alpha_df.pivot("alphaij", "alphaji", "I2")

I12_alpha_df = pd.DataFrame(dict(alphaij = a_ij, alphaji = a_ji, I12 = i_ij))
I12_to_alphas = I12_alpha_df.pivot("alphaij", "alphaji", "I12")

D_alpha_df = pd.DataFrame(dict(alphaij = a_ij, alphaji = a_ji, deceased = dead))
D_to_alphas = D_alpha_df.pivot("alphaij", "alphaji", "deceased")

########################	PLOT 	###########################################

############################	increase font size for plotting		#######################

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)  # pass in the font dict as kwargs

############## heatmap
ax = sns.heatmap(I1_to_alphas, vmin = -1, xticklabels=I1_to_alphas.columns.values.round(2),
                 yticklabels=I1_to_alphas.index.values.round(2), cbar_kws={'label': 'I1'})
plt.show()

############## heatmap
ax = sns.heatmap(I2_to_alphas, vmin=-1, xticklabels=I2_to_alphas.columns.values.round(2),
                 yticklabels=I2_to_alphas.index.values.round(2), cbar_kws={'label': 'I2'})
plt.show()

############## heatmap
ax = sns.heatmap(I12_to_alphas, vmin=-1, xticklabels=I12_to_alphas.columns.values.round(2),
                 yticklabels=I12_to_alphas.index.values.round(2), cbar_kws={'label': 'I12'})
plt.show()

############## heatmap
ax = sns.heatmap(D_to_alphas, xticklabels=D_to_alphas.columns.values.round(2),
                 yticklabels=D_to_alphas.index.values.round(2), cbar_kws={'label': 'deceased hosts'})
plt.show()
