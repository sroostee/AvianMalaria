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
c_delta1 = np.arange(0.0001,0.01,0.001)
c_delta2 = np.arange(0.0001,0.01,0.001)

s = []
i_1 = []
i_2 = []
i_ij = []
dead = []
c_d1 = []
c_d2 = []

for c_1 in c_delta1:
	for c_2 in c_delta2:
		c_d1.append(c_1)
		c_d2.append(c_2)
		#calculate copies of each strain in the host
		# number of copies of both strains
		n1 = inHost.K1
		n2 = inHost.K2
		#at equilibrium
		n1_12 = (inHost.K1 - inHost.alpha12 * inHost.K2)/(1- inHost.alpha12 *inHost.alpha21) 
		n2_12 = (inHost.K2 - inHost.alpha21 * inHost.K1)/(1- inHost.alpha12 *inHost.alpha21) 
		#calculate the number of hosts with a double infection
		host_out = odeint(host.eq_sys, host.y0, host.time, args =(c_1, c_2, 
			inHost.c_beta, inHost.h, n1, n2, n1_12, n2_12, host.mu, host.labda, host.l))
		S, I_1, I_2, I_12, dead_host = host_out.T
		S = S[-1]
		I_1 = I_1[-1]
		I_2 = I_2[-1]
		I_12 = I_12[-1]
		dead_host = dead_host[-1]
		s.append(S)
		i_1.append(I_1)
		i_2.append(I_2)
		i_ij.append(I_12)
		dead.append(dead_host)


S_cdelta_df = pd.DataFrame(dict(c_delta1 = c_d1, c_delta2 = c_d2, S = s))
S_to_cdeltas = S_cdelta_df.pivot("c_delta1", "c_delta2", "S")

I1_cdelta_df = pd.DataFrame(dict(c_delta1 = c_d1, c_delta2 = c_d2, I1 = i_1))
I1_to_cdeltas = I1_cdelta_df.pivot("c_delta1", "c_delta2", "I1")

I2_cdelta_df = pd.DataFrame(dict(c_delta1 = c_d1, c_delta2 = c_d2, I2 = i_2))
I2_to_cdeltas = I2_cdelta_df.pivot("c_delta1", "c_delta2", "I2")

I12_cdelta_df = pd.DataFrame(dict(c_delta1 = c_d1, c_delta2 = c_d2,  I12 = i_ij))
I12_to_cdeltas = I12_cdelta_df.pivot("c_delta1", "c_delta2", "I12")

D_cdelta_df = pd.DataFrame(dict(c_delta1 = c_d1, c_delta2 = c_d2,  deceased = dead))
D_to_cdeltas = D_cdelta_df.pivot("c_delta1", "c_delta2", "deceased")

########################	PLOT 	###########################################

############################	increase font size for plotting		#######################

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)  # pass in the font dict as kwargs

############## heatmap
ax = sns.heatmap(S_to_cdeltas, yticklabels=S_to_cdeltas.index.values.round(4), cbar_kws={'label': 'S'})
plt.show()

############## heatmap
ax = sns.heatmap(I1_to_cdeltas, yticklabels=I2_to_cdeltas.index.values.round(4), cbar_kws={'label': 'I1'})
plt.show()

############## heatmap
ax = sns.heatmap(I2_to_cdeltas, yticklabels=I2_to_cdeltas.index.values.round(4), cbar_kws={'label': 'I2'})
plt.show()

# ############## heatmap
ax = sns.heatmap(I12_to_cdeltas, xticklabels=I12_to_cdeltas.columns.values.round(4),
                 yticklabels=I12_to_cdeltas.index.values.round(4), cbar_kws={'label': 'I12'})
plt.show()

ax = sns.heatmap(D_to_cdeltas, xticklabels=D_to_cdeltas.columns.values.round(4),
                 yticklabels=D_to_cdeltas.index.values.round(4), cbar_kws={'label': 'deceased hosts'})
plt.show()





