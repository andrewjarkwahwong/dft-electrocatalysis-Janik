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
    
# Initial State Ex: Gas phase + Bare Surface
e_in = -221.416 #Energy of Initial State (eV)
dm_in = 0.03877 #Dipole moment of Initial state (eA)
polar_in_un = 4.549
polar_bare = 4.549
polar_in =polar_in_un - polar_bare


# Final State Ex: Adsorbed State
e_fin = -227.149199 #Energy of Final State (eV)
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
u_low = -2 #Lower Potential V-NHE
u_high = 1 #Upper Potential V-NHE


# EDL Model Parameters
er = 13 #Relative permittivity (Dielectric Constant)
e_vac = 0.00553 #Vacuum potential
d = 3 #Helmholtz EDL Width in Angstrom
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
e=er*e_vac #complex permittivity

u = np.linspace(u_low,u_high,25)
u_prime = u - u_pzc



#Faradaic


# Model 1A: Free Energy Change at U_PZC
g_1a = e_fin - e_in + g_solv + wf_bare - vac_nhe

# Model 1B: Potenital dependent Free Energy Change assuming Beta = 1 (Full e- transfer)
g_1b = g_1a+u_prime


#Model 2A: Incorporate Capacitance Charge
C = er*a*e_vac/d # Predicted Capacitiance by Helmholtz Model
diff_dm_sq = (dm_fin ** 2 - dm_in ** 2) 
# 0th order capacitance
C_0 = -0.5 * (diff_dm_sq) / (C * d ** 2)
#  1st order wrt U'
diff_dm = round(dm_fin - dm_in, 2)
C_const_1 = diff_dm / d
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



#Model 2C: Incorporating Polarizability (Induced Dipole-Field Terms)
diff_dm_polar_sq = (polar_fin*(dm_fin)*dm_fin-polar_in*dm_in*dm_in)
diff_a_dm = polar_fin*dm_fin-polar_in*dm_in
# Polarizability 0th order
p_0_d = 2*(e**2)*(a**2)*(d**2)
p_0 = diff_dm_polar_sq/p_0_d
# Polarizability 1st order
p_1_d=e*a*d**2
p_1 = -u_prime*diff_a_dm/p_1_d
# Polarizability 2nd order
p_2 = 0.5*(u_prime**2)*(diff_polar)*(1/d)*(1/d)
# Polarizability total
p_total = p_0 + p_1 + p_2 
#Calculate G_2C
g_2c=g_2b+p_total

    
# %% Free Energy Change vs Applied Potential across different EDL Models

#Calculate the roots
Model =['1b','2a','2b']
for i in Model:
    globals()['TL_'+i]=np.polyfit(u_prime,globals()['g_'+i],1) #Solve coefficents
    globals()['TL_'+i+'_1'] = np.poly1d(globals()['TL_'+i]) #Put coefficients into a polynomial
#2c model roots
TL_2c = np.polyfit(u_prime,g_2c,2)
TL_2c_1 = np.poly1d(TL_2c)



# Create Figure
fig= plt.figure(figsize=(12,10))
ap2=fig.add_subplot(111) #rows,columns,subgraph
ap2.plot(u_prime,TL_1b_1(u_prime),c='k',label="Model 1A:No EDL",linewidth=5)
ap2.plot(u_prime,TL_2a_1(u_prime),c='r',label="Model 2A:C",linewidth=5)
ap2.plot(u_prime,TL_2b_1(u_prime),c='b',label="Model 2B:C+$\mu$",linewidth=5)
ap2.plot(u_prime,TL_2c_1(u_prime),c='g',label=r"Model 2C:C+$\mu$+$\alpha$",linewidth=5)


#Figure Parameters
ap2.set_facecolor('white')
ap2.patch.set_edgecolor('black')  
ap2.patch.set_linewidth('5')  
ap2.grid(False)
plt.legend(loc='best',fontsize=28)
plt.xlabel('U-U$_{pzc}$ (V-NHE)', fontweight = 'bold', fontsize = 32)
plt.ylabel('Free Energy Change (eV)', fontweight = 'bold', fontsize = 32)
plt.title('Effects of EDL Models on $\Delta$G vs Applied Potential', fontweight = 'bold', fontsize = 36)

#Set the Axis Parameters
g_data=(g_1b,g_2a,g_2b,g_2c)
gmax=max(np.concatenate(g_data))
gmin=min(np.concatenate(g_data))

plt.axis([u_prime[0]-0.5,u_prime[-1]+0.5,gmin-0.5,gmax+0.5])
ap2.tick_params(axis='x', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
ap2.tick_params(axis='y', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)



#Text Parameters
plt.text(0.75 + u_prime[-1], gmax, "$\Delta\mu$ = " + str(round(diff_dm, 2)) + "eÅ",fontsize=28)
plt.text(0.75 + u_prime[-1], gmax-0.5, r"$\Delta$$\alpha$ = " + str(round(diff_polar, 2)) + "eÅ$^2$V$^{-1}$",fontsize=28)
plt.text(0.75+u_prime[-1],gmax-1,"$\epsilon_r$ = "+str(er),fontsize=28)
plt.text(0.75+u_prime[-1],gmax-1.5,"d$_{EDL}$ = " +str(d) + " Å",fontsize=28)
plt.text(0.75+u_prime[-1],gmax-2,"U$_{pzc}$ = " +str(round(u_pzc, 2)) + " V-NHE",fontsize=28)
plt.text(0.75+u_prime[-1],gmax-2.5,"$\Delta$G$_{1A}$=%.2fU'+%.2f"%(TL_1b[0],TL_1b[1]),c='k',fontsize=28)
plt.text(0.75+u_prime[-1],gmax-3,"$\Delta$G$_{2A}$=%.2fU'+%.2f"%(TL_2a[0],TL_2a[1]),c='r',fontsize=28)
plt.text(0.75+u_prime[-1],gmax-3.5,"$\Delta$G$_{2B}$=%.2fU'+%.2f"%(TL_2b[0],TL_2b[1]),c='b',fontsize=28)
plt.text(0.75+u_prime[-1],gmax-4,"$\Delta$G$_{2C}$=%.2fU'^2+%.2fU'+%.2f"%(TL_2c[0],TL_2c[1],TL_2c[2]),c='g',fontsize=28)

plt.show()

#Table 
print("Dielectric Constant =","{:.2f}".format(er), "EDL Width =","{:.2f}".format(d),"Å")
print("Dipole Moment Change =","{:.2f}".format(diff_dm),"eÅ", "Polarizability Change =","{:.2f}".format(diff_polar),"eÅ^2V^-1")
table_U0 = [['Model','Delta_G vs U-Upzc Equation'],['1B',"%.2fU'+%.2f"%(TL_1b[0],TL_1b[1])],['2A',"%.2fU'+%.2f"%(TL_2a[0],TL_2a[1])],['2B',"%.2fU'+%.2f"%(TL_2b[0],TL_2b[1])],['2C',"%.2fU'^2+%.2fU'+%.2f"%(TL_2c[0],TL_2c[1],TL_2c[2])]]
print(tb(table_U0,headers='firstrow',tablefmt = 'grid'))


# %% Contribution of each EDL component

# Create a figure with 3 subplots in a column
fig, axs = plt.subplots(3, 1, figsize=(8, 15), sharex=True)

# Plot Capacitance
axs[0].bar(u_prime, c_total, 0.05, label="Capacitance", color='red',zorder=0)

# Plot Dipole-Field
axs[1].bar(u_prime, dm_total, 0.05, label="Dipole-Field", color='blue',zorder=0)

# Plot Polarizability
axs[2].bar(u_prime, p_total, 0.05, label="Polarizability", color='green')
axs[2].set_xlabel("U-U$_{pzc}$ (V-NHE)", fontweight='bold', fontsize=32,zorder=0)

#Figure Parameters
for i in range(3):
    axs[i].set_facecolor('white')
    axs[i].patch.set_edgecolor('black')  
    axs[i].patch.set_linewidth('5')
    axs[i].set_ylabel("Free Energy Change (eV)", fontweight='bold', fontsize=20)
    axs[i].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
    axs[i].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
    axs[i].tick_params(axis='y',which='minor', length=5,width=4,direction='in')
    axs[i].tick_params(axis='x',which='minor', length=5,width=4,direction='in')
    axs[i].legend(loc='best', fontsize=24)
    axs[i].set_xlim([u_prime[0]-0.5,u_prime[-1]+0.5])
    axs[i].grid(False)

    total_edl=c_total+dm_total+p_total
    if max(total_edl) <= 0:
        axs[i].set_ylim(min(total_edl)-0.5,)
    else:
        axs[i].set_ylim(min(total_edl)-0.5,max(total_edl))



# Set common title
title_string = "Decompartmentalizing EDL Effects at $\epsilon_r$ = {:.2f} and d = {:.2f} Å".format(er, d)
fig.suptitle(title_string, fontweight='bold', fontsize=18)

# Adjust layout for better spacing
plt.tight_layout()

# Show the plot
plt.show()

#%% Total EDL COntribution

#Plot Total
total_edl=c_total+dm_total+p_total
fig= plt.figure(figsize=(12,10))
ap2=fig.add_subplot(111) #rows,columns,subgraph
plt.bar(u_prime,total_edl,0.1,color='blue')

colors = sns.color_palette('deep', n_colors=8)

#Figure Parameters
ap2.set_facecolor('white')
ap2.patch.set_edgecolor('black')  
ap2.patch.set_linewidth('5')  
ap2.grid(False)
plt.xlabel("U-U$_{pzc}$ (V-NHE)", fontweight = 'bold', fontsize = 32)
plt.ylabel("Free Energy Change(eV)", fontweight = 'bold', fontsize = 32)
plt.title("Total Contribution of EDL Effects" , fontweight = 'bold', fontsize = 32)


#Axis Set
ap2.tick_params(axis='x', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(autominorLocator())
ap2.tick_params(axis='y', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
ap2.tick_params(axis='y',which='minor', length=5,width=4,direction='in')
ap2.tick_params(axis='x',which='minor', length=5,width=4,direction='in')
total = np.abs(c_total) + np.abs(dm_total) + np.abs(p_total)
plt.xlim([u_prime[0]-0.5,u_prime[-1]+0.5])

if max(total_edl) <= 0:
    plt.ylim(min(total_edl)-0.5,0.5)
else:
    plt.ylim(min(total_edl)-0.5,max(total_edl)+0.5)
plt.show()



