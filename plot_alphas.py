#!/usr/bin/env python3

"""
date: 01/04/2019
author. S.J. Roostee

Visualisation of the effect of the varying competition coefficient values 
in Lotka-Volterra dynamics on the host population

"""
###########################		module import 		###############################

import numpy as np
from scipy.integrate import odeint
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import withinHost_model as inHost
import host_model as host

######################### generate alpha values ############################
alpha_ij = np.arange(0.1,1.0,0.1)
alpha_ji = np.arange(0.1,1.0,0.1)

s = []
i_1 = []
i_2 = []
i_ij = []
dead_inf = []
dead_normal = []
dead_total = []
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
		S, I_1, I_2, I_12, dead_i, dead_host = host_out.T
		S = S[-1]
		I_1 = I_1[-1]
		I_2 = I_2[-1]
		I_12 = I_12[-1]
		dead_i = dead_i[-1]
		dead_host = dead_host[-1]
		s.append(S)
		i_1.append(I_1)
		i_2.append(I_2)
		i_ij.append(I_12)
		dead_inf.append(dead_i)
		dead_normal.append(dead_host)
		dead_total.append(dead_host + dead_i)


S_alpha_df = pd.DataFrame(dict(alpha12 = a_ij, alpha21 = a_ji, S = s))
S_to_alphas = S_alpha_df.pivot("alpha12", "alpha21", "S")

I1_alpha_df = pd.DataFrame(dict(alpha12 = a_ij, alpha21 = a_ji, I1 = i_1))
I1_to_alphas = I1_alpha_df.pivot("alpha12", "alpha21", "I1")

I2_alpha_df = pd.DataFrame(dict(alpha12 = a_ij, alpha21 = a_ji, I2 = i_2))
I2_to_alphas = I2_alpha_df.pivot("alpha12", "alpha21", "I2")

I12_alpha_df = pd.DataFrame(dict(alpha12 = a_ij, alpha21 = a_ji, I12 = i_ij))
I12_to_alphas = I12_alpha_df.pivot("alpha12", "alpha21", "I12")

dead_norm_frac= [di/dt for di,dt in zip(dead_normal,dead_total)]

D_alpha_df = pd.DataFrame(dict(alpha12 = a_ij, alpha21 = a_ji, deceased = dead_total))
D_to_alphas = D_alpha_df.pivot("alpha12", "alpha21", "deceased")

dead_inf_fraction = [di/dt for di,dt in zip(dead_inf,dead_total)]

Df_alpha_df = pd.DataFrame(dict(alpha12 = a_ij, alpha21 = a_ji, deceased_fraction = dead_inf_fraction))
Df_to_alphas = Df_alpha_df.pivot("alpha12", "alpha21", "deceased_fraction")

########################	PLOT 	###########################################

############################	increase font size for plotting		#######################

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)  # pass in the font dict as kwargs

############## heatmap
ax = sns.heatmap(S_to_alphas, vmin = -1, xticklabels=S_to_alphas.columns.values.round(2),
                 yticklabels=S_to_alphas.index.values.round(2), cbar_kws={'label': 'S'})
plt.show()

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

############## heatmap
ax = sns.heatmap(Df_to_alphas, xticklabels=Df_to_alphas.columns.values.round(2),
                 yticklabels=Df_to_alphas.index.values.round(2), cbar_kws={'label': 'deceased hosts through infection'})
plt.show()
