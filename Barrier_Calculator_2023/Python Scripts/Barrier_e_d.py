# -*- coding: utf-8 -*-
"""
Created on Friday Feb 15 2024

Analytical GC framework for NH reduction on Rh (111) 

This code will allow you to do the following: 
-Investigate the compartmentalization of EDL effects at a given er and d
-Quantify the Sensitivity of er and d on potential-dependent activation energies
-Compartmentalize each EDL effect based on different er and d

The code uses pandas to read from excel (which I find easier to read data and can incorporate list of data easily)
You will need the aGCDFT_ed.xlsx excel sheet and this script to run this code. 
The excel formatting is important for pandas to read so maintain the same cell/array sizes

Feel free to reach out below if you have any questions.

@author: Andrew Jark-Wah Wong (Email: ajwongphd@gmail.com)
"""
# %% Inputs
import PySimpleGUI as sg
import ipywidgets 
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate as tb
import pandas as pd
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator
# Sheet name
sheet = 'NH_Rh' # TMA.py or TMA4x4.py
#Import Data
df = pd.read_excel("C:/Users/dreww/OneDrive/Documents/GitHub/dft-electrocatalysis-Janik/Barrier_Calculator_2023/Excel Sheets/aGCDFT_ed.xlsx", sheet_name=sheet)

M = ['NH* to NH2*'] #Data Name

#
a = df['area'].tolist() # Area of the cell (A^2)
u_pzc = df['upzc'].tolist() #Potential of zero charge of the bare metal surface. Here I use the workfunction of bare metal surface (V-NHE)

# Initial State
e_ref = df['G_In'].tolist() #Energy of Reference State (ex: adsorbate without H*) (eV)
e_H2 = df['G_H2'].tolist() #Energy of H2 gas (eV)
e_in = [i + y for i,y in zip(e_ref,e_H2)] # Energy of Initial State (eV)

dm_in = df['DM_In'].tolist()  #Dipole moment of Initial state (eA)
polar_in_un = df['Polar_In'].tolist() #Polarizability of the reference state (eA^2V-1)
polar_bare = df['Polar_Bare'].tolist() #Polarizability of the bare surface (eA^2V-1)
polar_in =[i - y for i,y in zip(polar_in_un,polar_bare)] #Corrected polarizability of the reference state adsorbated (eA^2V-1)

# Final State Ex: Adsorbed State
e_fin = df['G_Fin'].tolist() #Energy of Initial State (eV)
dm_fin = df['DM_Fin'].tolist()  #Dipole moment of Initial state (eA)
polar_fin_un = df['Polar_Fin'].tolist()  #Polarizability of the Final state (eA^2V-1)
polar_fin = [i - y for i,y in zip(polar_fin_un,polar_bare)] #Corrected polarizability of the Final statestate adsorbated (eA^2V-1)

#Potential Range of Interest
vac_nhe = 4.6 #Converting from V-Abs to V-NHE typically varies from 4.2 to 5 V
u_low = -2.5 #Lower Potential V-NHE
u_high = 1 #Upper Potential V-NHE
diff_dm_sq = [i ** 2 - y ** 2 for i,y in zip(dm_fin,dm_in)] 


# EDL Model Parameters
e_vac = 0.00553 #Vacuum potential in eV^-1A^-1
g_solv = df['G_Solv'].tolist() #solvation free energy change of the reaction (difference of solvation free energy of products - reactants)

#Single Dielectric constants for Figure 1
er = [2] #Relative permittivity (Dielectric Constant)
d = [3] #Helmholtz EDL Width in Angstrom

#Set of dielectric constants and d's for Figure 2 and 3
#er = [1,2,4,8,13,78.4] #Relative permittivity (Dielectric Constant)
#d = [3,4.5,6,100] #Helmholtz EDL Width in Angstrom



# %% The Math

#Dipole moment and Polarizability Changes
diff_dm = [i - y for i,y in zip(dm_fin,dm_in)]
diff_polar = [i - y for i,y in zip(polar_fin,polar_in)]
diff_dm_polar_sq = [y * x * x - z * a * a for x,y,z,a in zip(dm_fin,polar_fin,dm_in,polar_in)]
diff_a_dm = [x * y - z* a for x,y,z,a in zip(dm_fin,polar_fin,dm_in,polar_in)]
             
#Potential 
u = np.linspace(u_low,u_high,25)
u_prime = u - u_pzc

# Define lists of values for er and d
er_values = er  # Relative permittivity (Dielectric Constant)
d_values = d  # Helmholtz EDL Width in Angstrom

# Initialize results dictionary
resultsEDL = {}
# Nested loops to iterate over er and d values
for er_val in er_values:
    for d_val in d_values:
        for i in range(len(M)):

            # Calculate g_1a for the current iteration
            g_1a = e_fin[i] - e_in[i] + g_solv[i] + u_pzc[i]
            
            # Calculate g_1b for the current iteration
            g_1b = g_1a + u_prime
            
            # Calculate capacitance terms
            C = er_val * a[i] * e_vac / d_val
            diff_dm_sq = (dm_fin[i] ** 2 - dm_in[i] ** 2) 
            C_0 = -0.5 * (diff_dm_sq) / (C * d_val ** 2)
            C_const_1 = diff_dm[i] / d_val
            C_1 = C_const_1 * u_prime
            c_total = C_0 + C_1
            
            # Calculate g_2a for the current iteration
            g_2a = g_1b + c_total
            
            # Calculate dipole-field terms
            dm_0 = 2 * C_0
            dm_1 = u_prime * C_const_1
            dm_total = dm_0 + dm_1
            g_2b = g_2a + dm_total
            
                        
            # Model 2C: Incorporating Polarizability (Induced Dipole-Field Terms)
            diff_dm_polar_sq = (polar_fin[i] * dm_fin[i] * dm_fin[i] - polar_in[i] * dm_in[i] * dm_in[i])
            diff_a_dm = polar_fin[i] * dm_fin[i] - polar_in[i] * dm_in[i]
            # Polarizability 0th order
            p_0_d = 2 * ((er_val*e_vac) ** 2) * (a[i] ** 2) * (d_val** 2)
            p_0 = diff_dm_polar_sq / p_0_d
            # Polarizability 1st order
            p_1_d = er_val*e_vac * a[i] * d_val ** 2
            p_1 = -1*u_prime * diff_a_dm / p_1_d
            # Polarizability 2nd order
            p_2 = 0.5 * (u_prime ** 2) * (diff_polar[i]) * (1 / d_val) * (1 / d_val)
            # Polarizability total
            p_total = p_0 + p_1 + p_2
            # Calculate G_2C
            g_2c = g_2b + p_total
            EDL_total = p_total + dm_total+c_total

            # Store the result in the dictionary
            resultsEDL[(er_val, d_val)] = {'u': u, 'g_1b':g_1b,'g_2a':g_2a,'g_2b':g_2b, 'g_2c': g_2c,'far':u_prime,'cap':c_total,'dm':dm_total,'p':p_total}

# %% Models based on single er and d

#If you specify multiple dielectric constant, they will ALL be plotted so specify here either one of the stored dictionary data 
#Or specify a single dielectric constant and d

import itertools           


palette = ['Blues', 'flare', 'crest']
colors = sns.color_palette(palette[1], n_colors=6)
color_cycle = itertools.cycle(colors)
# Create a single subplot


fig, ax = plt.subplots(figsize=(16,12))
 
# Plot u against g_2c in the current subplot
for (er_val, d_val), data in resultsEDL.items():
    ax.plot(data['u'], data['g_1b'],label="Faradaic No EDL",linewidth=5, alpha=1, color='black')
    ax.plot(data['u'], data['g_2a'],label=r"C",linewidth=5, alpha=1, color='red')
    ax.plot(data['u'], data['g_2b'],label=r"C+$\mu$",linewidth=5, alpha=1, color='blue' )
    ax.plot(data['u'], data['g_2c'],label=r"C+$\mu$+$\alpha$",linewidth=5, alpha=1, color='green')
    
    
    #Calculate the roots
    Model =['1b','2a','2b']
    for i in Model:
        globals()['TL_'+i]=np.polyfit(u,data['g_'+i],1) #Solve coefficents
        globals()['TL_'+i+'_1'] = np.poly1d(globals()['TL_'+i]) #Put coefficients into a polynomial
    #2c model roots
    TL_2c = np.polyfit(data['u'],data['g_2c'],2)
    TL_2c_1 = np.poly1d(TL_2c)

gmax = max(data['g_1b'])
gmin = min(data['g_1b'])

plt.text(0.1+u[-1],gmax-.9,"$\epsilon_r$ = "+str(er[0]),fontsize=32)
plt.text(0.1+u[-1],gmax-1.2,"d$_{EDL}$ = " +str(d[0]) + r" $\AA$",fontsize=32)
plt.text(0.1+u[-1],gmax-1.5,r"$\Delta$G$_{1B}$=%.2fU+%.2f"%(TL_1b[0],TL_1b[1]),c='k',fontsize=28)
plt.text(0.1+u[-1],gmax-1.8,r"$\Delta$G$_{2A}$=%.2fU+%.2f"%(TL_2a[0],TL_2a[1]),c='r',fontsize=28)
plt.text(0.1+u[-1],gmax-2.1,r"$\Delta$G$_{2B}$=%.2fU+%.2f"%(TL_2b[0],TL_2b[1]),c='b',fontsize=28)
plt.text(0.1+u[-1],gmax-2.4,r"$\Delta$G$_{2C}$=%.2fU$^2$+%.2fU+%.2f"%(TL_2c[0],TL_2c[1],TL_2c[2]),c='g',fontsize=28)


# Set labels for the subplot
ax.set_xlabel("U (V-SHE)", fontsize=37)
ax.set_ylabel("Activation Energy (V-SHE)", fontsize=36)

# Set limits and ticks for the subplot
ax.set_xlim([u[0],u[-1]])
ax.set_ylim([gmin, gmax])
#ax.xaxis.set_ticks(np.arange(-1.5, 0.25, 0.25))

# Set axis parameters for the subplot
ax.tick_params(axis='both', labelsize=28, width=4, colors='black', direction="in", grid_color='black', which='major',
               length=10, pad=15)
ax.tick_params(axis='both', which='minor', length=6, width=5, direction='in')

# Set ticks on the top and right sides
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

# Set spines visible on all sides
ax.spines['top'].set_linewidth(4)
ax.spines['right'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.spines['left'].set_linewidth(4)

ax.tick_params(axis='y', which='minor', length=6, width=5, direction='in')
ax.tick_params(axis='x', which='minor', length=6, width=5, direction='in')

ax.set_facecolor('white')
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)
minor_locator1 = AutoMinorLocator(1)
ax.yaxis.set_minor_locator(minor_locator1)
ax.patch.set_edgecolor('black')
ax.patch.set_linewidth(5)

#plt.tight_layout()
plt.grid(False)
plt.legend(loc='best',fontsize=32,ncol=1)
plt.show()

# %% Sensitivity to Barrier based on er and d
# Change the er and d

import itertools           
palette = ['Blues', 'flare', 'crest']
colors = sns.color_palette(palette[1], n_colors=6)
color_cycle = itertools.cycle(colors)
# Create a single subplot

fig, ax = plt.subplots(figsize=(16,12))


# Enumerate over combinations of er_val and d_val
for (er_val, d_val), data in resultsEDL.items():
    if d_val == 3:
        linestyle = '-'
    elif d_val == 4.5:
        linestyle = '-.'
    elif d_val == 6:
        linestyle = ':'
    else:
        linestyle = ':'

    # Get the next color from the color cycle
    color = next(color_cycle)

    # Plot u against g_2c in the current subplot
    line, = ax.plot(data['u'], data['g_2c'],
                    label=fr"e$_r$={er_val}, d={d_val} $\AA$",
                    linewidth=5, alpha=1, color=color, linestyle=linestyle)
    
    



# Set labels for the subplot
ax.set_xlabel("U (V-SHE)", fontweight='bold', fontsize=32)
ax.set_ylabel("Adsorption Free Energy (V-SHE)", fontweight='bold', fontsize=30)

# Set limits and ticks for the subplot
ax.set_xlim([u[0], u[-1]])
ax.set_ylim([gmin, gmax])
#ax.xaxis.set_ticks(np.arange(-1.5, 0.25, 0.25))

# Set axis parameters for the subplot
ax.tick_params(axis='both', labelsize=28, width=4, colors='black', direction="in", grid_color='black', which='major',
               length=10, pad=15)
ax.tick_params(axis='both', which='minor', length=6, width=5, direction='in')

# Set ticks on the top and right sides
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

# Set spines visible on all sides
ax.spines['top'].set_linewidth(4)
ax.spines['right'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.spines['left'].set_linewidth(4)

ax.tick_params(axis='y', which='minor', length=6, width=5, direction='in')
ax.tick_params(axis='x', which='minor', length=6, width=5, direction='in')

ax.set_facecolor('white')
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)
minor_locator1 = AutoMinorLocator(1)
ax.yaxis.set_minor_locator(minor_locator1)
ax.patch.set_edgecolor('black')
ax.patch.set_linewidth(5)

plt.tight_layout()
plt.grid(False)
plt.legend(loc='best',fontsize=18,ncol=2)
plt.show()


# %% Decompartmentalization based on sets of er and d
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator
import itertools

palette = ['Blues', 'flare', 'crest']
colors = sns.color_palette(palette[1], n_colors=6)
color_cycle = itertools.cycle(colors)

# Create a single subplot
fig, axs = plt.subplots(1, 4, figsize=(24, 8),sharey=False)

# Enumerate over combinations of er_val and d_val
for (er_val, d_val), data in resultsEDL.items():
    if d_val == 3:
        linestyle = '-'
    elif d_val == 4.5:
        linestyle = '--'
    elif d_val == 6:
        linestyle = '-.'
    else:
        linestyle = ':'

    # Get the next color from the color cycle
    color = next(color_cycle)

    # Plot u against g_2c in the current subplot
    axs[0].plot(data['u'], data['far'],
                label=fr"e$_r$={er_val}, d={d_val}",
                linewidth=5, alpha=1, color=color, linestyle=linestyle)
    axs[1].plot(data['u'], data['cap'],
                label=fr"e$_r$={er_val}, d={d_val} $\AA$",
                linewidth=5, alpha=1, color=color, linestyle=linestyle)

    axs[2].plot(data['u'], data['dm'],
                label=fr"e$_r$={er_val}, d={d_val} + $\AA$",
                linewidth=5, alpha=1, color=color, linestyle=linestyle)

    axs[3].plot(data['u'], data['p'],
                label=fr"e$_r$={er_val}, d={d_val}",
                linewidth=5, alpha=1, color=color, linestyle=linestyle)

# Set labels for the subplots
axs[0].set_xlabel("U (V-SHE)", fontweight='bold', fontsize=28)
axs[0].set_ylabel("Faradaic (eV)", fontweight='bold', fontsize=32)
axs[1].set_xlabel("U (V-SHE)", fontweight='bold', fontsize=28)
axs[1].set_ylabel("Capacitive (eV)", fontweight='bold', fontsize=32)
axs[2].set_xlabel("U (V-SHE)", fontweight='bold', fontsize=28)
axs[2].set_ylabel("Dipole-Field (eV)", fontweight='bold', fontsize=32)
axs[3].set_xlabel("U (V-SHE)", fontweight='bold', fontsize=28)
axs[3].set_ylabel("Induced Dipole-Field (eV)", fontweight='bold', fontsize=28)

# Set limits and ticks for the subplots
for ax in axs:
    ax.set_xlim([u[0], u[-1]])
    ax.set_ylim([gmin, gmax])
    ax.tick_params(axis='both', labelsize=16, width=2, colors='black', direction="in", grid_color='black', which='major',
               length=8, pad=10)
    ax.tick_params(axis='both', which='minor', length=4, width=3, direction='in')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.set_facecolor('white')
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    minor_locator1 = AutoMinorLocator(1)
    ax.yaxis.set_minor_locator(minor_locator1)
    ax.patch.set_edgecolor('black')
    ax.patch.set_linewidth(5)

# Add legend to the first subplot
axs[1].legend(loc='best', fontsize=18, ncol=1)

plt.tight_layout()
plt.grid(False)
plt.show()
