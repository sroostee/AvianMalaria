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

######################### generate K values ############################
K_1 = np.arange(3,15,.5)
K_2 = np.arange(3,15,.5)

i_1 = []
i_2 = []
i_12 = []
dead = []
K_1_ij = []
K_2_ij = []

for k1 in K_1:
	for k2 in K_2:
		#store K values for dataframe
		K_1_ij.append(k1)
		K_2_ij.append(k2)
		#calculate copies of each strain in the host
		# number of copies of both strains
		n1 = k1
		n2 = k2
		#at equilibrium
		n1_12 = (k1 - inHost.alpha12 * k2)/(1- inHost.alpha12*inHost.alpha21) 
		n2_12 = (k2 - inHost.alpha21 * k1)/(1- inHost.alpha12*inHost.alpha21) 
		#calculate the number of hosts with a double infection
		host_out = odeint(host.eq_sys, host.y0, host.time, args =(inHost.c_delta1, inHost.c_delta2, 
			inHost.c_beta, inHost.h, n1, n2, n1_12, n2_12, host.mu, host.labda, host.l))
		S, I_1, I_2, I_12, dead_host = host_out.T
		I_1 = I_1[-1]
		I_2 = I_2[-1]
		I_12 = I_12[-1]
		dead_host = dead_host[-1]
		i_1.append(I_1)
		i_2.append(I_2)
		i_12.append(I_12)
		dead.append(dead_host)


I1_K_df = pd.DataFrame(dict(K1 = K_1_ij, K2 = K_2_ij, I12 = i_1))
I1_to_Ks = I1_K_df.pivot("K1", "K2", "I12")

I2_K_df = pd.DataFrame(dict(K1 = K_1_ij, K2 = K_2_ij, I12 = i_2))
I2_to_Ks = I2_K_df.pivot("K1", "K2", "I12")

I12_K_df = pd.DataFrame(dict(K1 = K_1_ij, K2 = K_2_ij, I12 = i_12))
I12_to_Ks = I12_K_df.pivot("K1", "K2", "I12")

D_K_df = pd.DataFrame(dict(K1 = K_1_ij, K2 = K_2_ij, deceased = dead))
D_to_Ks = D_K_df.pivot("K1", "K2", "deceased")

########################	PLOT 	###########################################

############################	increase font size for plotting		#######################

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)  # pass in the font dict as kwargs

############## heatmap
ax = sns.heatmap(data = I1_to_Ks, cbar_kws={'label': 'I1'})
plt.show()

ax = sns.heatmap(data = I2_to_Ks, cbar_kws={'label': 'I2'})
plt.show()

ax = sns.heatmap(data = I12_to_Ks, cbar_kws={'label': 'I12'})
plt.show()

ax = sns.heatmap(data = D_to_Ks, cbar_kws={'label': 'deceased hosts'})
plt.show()


