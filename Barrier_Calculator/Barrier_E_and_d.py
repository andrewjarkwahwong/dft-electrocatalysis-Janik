# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 00:33:07 2023

Compartmentalized

@author: Andrew Jark-Wah Wong (Email: ajwongphd@gmail.com)
"""

import PySimpleGUI as sg
import ipywidgets 
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate as tb
import pandas as pd
import seaborn as sns

# Initial State Ex: Gas phase + Bare Surface
e_in = -221.416 #Energy of Initial State (eV)
dm_in = 0.03877 #Dipole moment of Initial state (eA)
polar_in_un = 4.549
polar_bare = 4.549
polar_in =polar_in_un - polar_bare


# Final State Ex: Adsorbed State
e_fin =  -226.931535 #Energy of Final State (eV)
dm_fin = -1.518251 #Dipole Moment of Final State (eA)
polar_fin_un = 6.2654
polar_fin = polar_fin_un-polar_bare


# Cell Parameters
v = 1584.35 #Volume of Slab in (Angstrom^3)
h = 64.87335 #Height of Slab (Angstrom)
e_fermi = 1.599 #Fermi Energy (eV)
u_vac = 6.661  #Vacuum Potential (V-Abs)


#Potential Range of Interest
vac_nhe = 4.6 #Converting from V-Abs to V-NHE
u_low = -1.5 #Lower Potential V-NHE
u_high = 1 #Upper Potential V-NHE


# EDL Model Parameters
er = [2,3,4,6,8,13,78] #Relative permittivity (Dielectric Constant)
e_vac = 0.00553 #Vacuum potential
d = [3.5,4,5,6] #Helmholtz EDL Width in Angstrom
g_solv = 0 #solvation free energy change of the reaction

# %% Math

#Dipole moment and Polarizability Changes
diff_dm = dm_fin - dm_in
diff_polar = polar_fin - polar_in

#Constants for the cell
a=v/h # area of the slab (A^2)
wf_bare = u_vac-e_fermi #Work function of the bare metal surface
u_pzc = wf_bare - vac_nhe #Potential of zero charge of the bare surface

#Constants related to EDL model
e_vac = 0.00553 #vacuum permittivity
e= [i*e_vac for i in er] #complex permittivity

u = np.linspace(u_low,u_high,25)
u_prime = u - u_pzc



#Faradaic

answer = input(" Is this a Faradaic Reaction? Type YES or NO")
# Results dictionary to store the calculated values
results = {}

if answer.lower() == "yes":
    for er_val in er:
        for d_val in d:
            # Model 1A: Free Energy Change at U_PZC
            g_1a = e_fin - e_in + g_solv + wf_bare - vac_nhe
            
            # Model 1B: Potenital dependent Free Energy Change assuming Beta = 1 (Full e- transfer)
            g_1b = g_1a+u_prime
            
            
            #Model 2A: Incorporate Capacitance Charge
            C = er_val*a*e_vac/d_val # Predicted Capacitiance by Helmholtz Model
            diff_dm_sq = (dm_fin ** 2 - dm_in ** 2) 
            # 0th order capacitance
            C_0 = -0.5 * (diff_dm_sq) / (C * d_val ** 2)
            #  1st order wrt U'
            diff_dm = round(dm_fin - dm_in, 2)
            C_const_1 = diff_dm / d_val
            C_1 = C_const_1 * u_prime
            #  Total Capacitance
            c_total = C_0 + C_1 #  Total Capacitance
            # Free Energy Change to correct Model 1B with Charging
            g_2a = g_1b + c_total
            
            
            #Model 2B: Incorporating Dipole-Field Terms
            dm_0 = 2 * C_0
            dm_1 = u_prime * C_const_1
            dm_total = dm_0 + dm_1
            # Calculate G_2B to correct Model 2A with dipole-field terms
            g_2b = g_2a + dm_total
            
        
            # Model 2C: Incorporating Polarizability (Induced Dipole-Field Terms)
            diff_dm_polar_sq = (polar_fin * (dm_fin) * dm_fin - polar_in * dm_in * dm_in)
            diff_a_dm = polar_fin * dm_fin - polar_in * dm_in
            # Polarizability 0th order
            p_0_d = 2 * (er_val*e_vac ** 2) * (a ** 2) * (d_val ** 2)
            p_0 = diff_dm_polar_sq / p_0_d
            # Polarizability 1st order
            p_1_d = er_val*e_vac * a * d_val ** 2
            p_1 = -u_prime * diff_a_dm / p_1_d
            # Polarizability 2nd order
            p_2 = 0.5 * (u_prime ** 2) * (diff_polar) * (1 / d_val) * (1 / d_val)
            # Polarizability total
            p_total = p_0 + p_1 + p_2
            # Calculate G_2C
            g_2c = g_2b + p_total
            # Now you can use g_2c for further analysis or store the results as needed
            # Store the result in the dictionary
            results[(er_val, d_val)] = {'u_prime': u_prime, 'g_2c': g_2c}
            
# %% Plot 
# Plotting
colors = sns.color_palette('deep', n_colors=30)
fig = plt.figure(figsize=(18, 14))
ap2=fig.add_subplot(111) #rows,columns,subgraph

for i, ((er_val, d_val), data) in enumerate(results.items()):
    if d_val == 3:
        plt.plot(data['u_prime'], data['g_2c'], label=f'$\epsilon_r$={er_val}, d={d_val}', linewidth=5, color=colors[i+3])
    elif d_val ==3.5:
        plt.plot(data['u_prime'], data['g_2c'], label=f'$\epsilon_r$={er_val}, d={d_val}', linewidth=7.5, linestyle=':', color=colors[i+2])
    elif d_val ==4:
        plt.plot(data['u_prime'], data['g_2c'], label=f'$\epsilon_r$={er_val}, d={d_val}', linewidth=5, linestyle='dashdot', color=colors[i+1])

plt.xlabel("U-U$_{pzc}$ (V-NHE)", fontweight='bold', fontsize=56)
plt.ylabel("$G_{ads}$ (eV)", fontweight='bold', fontsize=56)
#plt.title("TMA Adsorption Free Energy Profile at different EDL Properties", fontweight='bold', fontsize=56)
plt.legend(fontsize=32,loc=(1.01, 0.1), ncol=1)
plt.xlim([u_prime[0]+0.45,0.])
plt.ylim([-5.8,-4.6])

#Axis Set
ap2.tick_params(axis='x', labelsize=42,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=25)
#ax.yaxis.set_minor_locator(autominorLocator())
ap2.tick_params(axis='y', labelsize=42,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=25)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
ap2.tick_params(axis='y',which='minor', length=5,width=4,direction='in')
ap2.tick_params(axis='x',which='minor', length=5,width=4,direction='in')
ap2.set_facecolor('white')
ap2.patch.set_edgecolor('black')  
ap2.patch.set_linewidth('5')  

plt.tight_layout()
plt.grid(False)
plt.show()