# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:24:03 2023
Beta Calculator 
Here the python script elucidate the sensitivity of the symmetry factor based on the four following parameters:
1. EDL Width (d)
2. Dielectric Constant
3. Magnitude of the Dipole Moment Change
4. Magnitude of the Polarizability Change     

Relevant Paper: 
1. (JACS Paper)
2. https://doi.org/10.1016/j.jcat.2024.115360

@author: Andrew Jark-Wah Wong
Email: ajwongphd@gmail.com
"""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.interpolate import interp1d
plt.style.use('seaborn')
colors = sns.color_palette('deep', n_colors=8)

# %% Dielectric Constant Varies

#Constants of the Cell and EDL
e_r = [1,2,8,13,78.4] #Relative Permittivity Change this
d = np.linspace(1,20,50)
e_vac = 0.00553 # Vacuum Permittivity 
u_pzc = 0.462 #Au 111 Upzc (V-SHE)
A = 64.87335 #Au 111

#Dipole Moment Change (eA)
dm_bare = 0.0 #Dipole Moment of the bare surface (eA)
dm_fin=[1] #Dipole Moment of the Final State (eA)

#Polarizability (eA^2V^-1)
polar_bare = 4.5 #Polarizability of the bare metal surface
polar_in = 4.5 - polar_bare #Polarazibility corrected of the adsorbate (none)
polar_fin_un = [4.5,4.75,5,5.25,6] #Polarizability of Final State
polar_fin = [i - polar_bare for i in polar_fin_un] #Polarizability

u=np.linspace(-1,1,50)

#Changes in Polarizability and Dipole Moments 
polar = np.array([i - polar_in for i in polar_fin])
dm = np.array([i - dm_bare for i in dm_fin])
a_dm= np.array([dm_fin[0] * y - dm_bare*polar_in for y in polar])
u_prime =np.array( [i-u_pzc for i in u])

#Beta varies with dm
beta_avg_dm = []  # Create a list to store separate lists of beta_avg

desired_index_dm = 0 # index of dm
desired_index_a = 1 # index of polar and a_dm

#Beta is calculated based on non-faradaic step 
for index in range(len(polar)):
    beta_avg = []  # Create a new list for each dm value
    for i in d: #polar and a_dm
        a_beta = polar[desired_index_a] / (i ** 2)
        b_beta = 2 * dm[desired_index_dm] / i - (a_dm[desired_index_a] / (e_r[index] * e_vac* A * i ** 2))
        beta = b_beta + a_beta * u_prime
        beta_avg_i = sum(beta) / len(beta)
        beta_avg.append(beta_avg_i)
    beta_avg_dm.append(beta_avg)  # Append the beta_avg list to beta_avg_lists

fig,ax = plt.subplots(figsize=(12,10))
#ax2=ax.twinx()
for i in range(len(polar)):
    ax.plot(d,beta_avg_dm[i],color=colors[i],label=r'$\epsilon$$_r$ = ' + str(e_r[i]),linewidth=5)
plt.axhline(y=0.5,color='red',linewidth=4,linestyle='--',label='\u03B2 = 0.5')
ax.set_facecolor('white')
ax.patch.set_edgecolor('black')  
ax.patch.set_linewidth(5)  

ax.text(3,0.925,r'$\Delta$$\alpha$ = ' + str(polar[desired_index_a]) + r'e$\AA$$^2$V$^-$',fontsize=28)
ax.text(3,0.85,r'$\Delta$$\mu$ = ' + str(dm[desired_index_dm])+r'e$\AA$',fontsize=28)


#labels
ax.set_ylabel('Symmetry Factor',fontsize=32,weight='bold',color='black')
ax.set_xlabel('Helmholtz EDL Width ($\AA$)',fontsize=32,weight='bold',color='black')
ax.grid(False)
ax.set_ylim(0,1)
ax.set_xlim(1,10)
ax.tick_params(axis='x', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(autominorLocator())
ax.tick_params(axis='y', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y',which='minor', length=5,width=4,direction='in')
ax.tick_params(axis='x',which='minor', length=5,width=4,direction='in')
ax.xaxis.grid(False, which='minor')
#ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_ticks(np.arange(1, 11, 1))
#ax.text(0.,0.65,'a)',fontsize=32,fontweight='bold')
legend=plt.legend(loc='best',fontsize=24,handlelength=2, handletextpad=1)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('black')
for text in legend.get_texts():
    text.set_color('black')
    text.set_weight('bold')
plt.show()

# %% DM Changes Sensitivity 

#Here the dipole moment changes are varied by varying the dipole moment of the final state (dm_fin). e_r and polar_fin need to be adjusted

#Constants of the Cell and EDL
e_r = [1] #Relative Permittivity Change this
d = np.linspace(1,20,50)
e_vac = 0.00553 # Vacuum Permittivity 
u_pzc = 0.462 #Au 111 Upzc (V-SHE)
A = 64.87335 #Au 111

#Dipole Moment Change
dm_bare = 0.0 #Dipole Moment of the bare surface (eA)
dm_fin=[0.0,0.25,0.5,0.75,1] #Dipole Moment of the Final State (eA)

#Polarizability
polar_bare = 4.5 #polarizability of the bare metal surface
polar_in = 4.5 - polar_bare #Polarazibility corrected of the adsorbate (none)
polar_fin_un = [4.5,4.75,5,5.25,6] #Polarizability of metal surface
polar_fin = [i - polar_bare for i in polar_fin_un] #Polarizability

u=np.linspace(-2,1,50)

#Changes in Polarizability and Dipole Moments 
polar = np.array([i - polar_in for i in polar_fin])
dm = np.array([i - dm_bare for i in dm_fin])
a_dm= np.array([x * y - dm_bare*polar_in for x,y in zip(polar_fin,dm_fin)])
u_prime =np.array( [i-u_pzc for i in u])

#Beta varies with dm
beta_avg_dm = []  # Create a list to store separate lists of beta_avg

#Beta is calculated based on non-faradaic step 
for dm_index in range(len(dm)):
    beta_avg = []  # Create a new list for each dm value
    for i in d:
        a_beta = polar[0] / (i ** 2)
        b_beta = 2 * dm[dm_index] / i - (a_dm[0] / (e_r[0] * e_vac* A * i ** 2))
        beta = b_beta + a_beta * u_prime
        beta_avg_i = sum(beta) / len(beta)
        beta_avg.append(beta_avg_i)
    beta_avg_dm.append(beta_avg)  # Append the beta_avg list to beta_avg_lists

fig,ax = plt.subplots(figsize=(12,10))
#ax2=ax.twinx()
for i in range(len(dm)):
    ax.plot(d,beta_avg_dm[i],color=colors[i],label = '$\Delta$$\mu$ = ' + str(dm[i]),linewidth=5)
plt.axhline(y=0.5,color='red',linewidth=4,linestyle='--',label='\u03B2 = 0.5')
ax.set_facecolor('white')
ax.patch.set_edgecolor('black')  
ax.patch.set_linewidth(5)  


ax.text(4.5,0.925,r'$\epsilon$$_r$ = ' + str(e_r[0]),fontsize=28)
ax.text(4.5,0.85,r'$\Delta$$\alpha$ = ' + str(polar[0])+r'e$\AA$$^2$V$^-$',fontsize=28)


#labels
ax.set_ylabel('Symmetry Factor',fontsize=32,weight='bold',color='black')
ax.set_xlabel('Helmholtz EDL Width ($\AA$)',fontsize=32,weight='bold',color='black')
ax.grid(False)
ax.set_ylim(0,1)
ax.set_xlim(1,10)
ax.tick_params(axis='x', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(autominorLocator())
ax.tick_params(axis='y', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y',which='minor', length=5,width=4,direction='in')
ax.tick_params(axis='x',which='minor', length=5,width=4,direction='in')
ax.xaxis.grid(False, which='minor')
#ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_ticks(np.arange(1, 11, 1))
#ax.text(0.,0.65,'a)',fontsize=32,fontweight='bold')
legend=plt.legend(loc='best',fontsize=24,handlelength=2, handletextpad=1)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('black')
for text in legend.get_texts():
    text.set_color('black')
    text.set_weight('bold')
plt.show()

# %% Polarizability Sensitivity 

#Constants of the Cell and EDL
e_r = [1] #Relative Permittivity Change this
d = np.linspace(1,20,50)
e_vac = 0.00553 # Vacuum Permittivity 
u_pzc = 0.462 #Au 111 Upzc (V-SHE)
A = 64.87335 #Au 111

#Dipole Moment Change
dm_bare = 0.0 #Dipole Moment of the bare surface (eA)
dm_fin=[0.25] #Dipole Moment of the Final State (eA)

#Polarizability
polar_bare = 4.5 #polarizability of the bare metal surface
polar_in = 4.5 - polar_bare #Polarazibility corrected of the adsorbate (none)
polar_fin_un = [4.5,4.75,5,5.25,6] #Polarizability of metal surface
polar_fin = [i - polar_bare for i in polar_fin_un] #Polarizability

u=np.linspace(-2,1,50)

#Changes in Polarizability and Dipole Moments 
polar = np.array([i - polar_in for i in polar_fin])
dm = np.array([i - dm_bare for i in dm_fin])
a_dm= np.array([dm_fin[0] * y - dm_bare*polar_in for y in polar])
u_prime =np.array( [i-u_pzc for i in u])

#Beta varies with dm
beta_avg_dm = []  # Create a list to store separate lists of beta_avg

#Beta is calculated based on non-faradaic step 
for dm_index in range(len(polar)):
    beta_avg = []  # Create a new list for each dm value
    for i in d:
        a_beta = polar[dm_index] / (i ** 2)
        b_beta = 2 * dm[0] / i - (a_dm[dm_index] / (e_r[0] * e_vac* A * i ** 2))
        beta = b_beta + a_beta * u_prime
        beta_avg_i = sum(beta) / len(beta)
        beta_avg.append(beta_avg_i)
    beta_avg_dm.append(beta_avg)  # Append the beta_avg list to beta_avg_lists

fig,ax = plt.subplots(figsize=(12,10))
#ax2=ax.twinx()
for i in range(len(polar)):
    ax.plot(d,beta_avg_dm[i],color=colors[i],label=r'$\Delta$$\alpha$ = ' + str(polar[i]) + r'e$\AA$$^2$V$^-$',linewidth=5)
plt.axhline(y=0.5,color='red',linewidth=4,linestyle='--',label='\u03B2 = 0.5')
ax.set_facecolor('white')
ax.patch.set_edgecolor('black')  
ax.patch.set_linewidth(5)  

ax.text(3,0.925,r'$\epsilon$$_r$ = ' + str(e_r[0]),fontsize=28)
ax.text(3,0.85,r'$\Delta$$\mu$ = ' + str(dm[0])+r'e$\AA$',fontsize=28)


#labels
ax.set_ylabel('Symmetry Factor',fontsize=32,weight='bold',color='black')
ax.set_xlabel('Helmholtz EDL Width ($\AA$)',fontsize=32,weight='bold',color='black')
ax.grid(False)
ax.set_ylim(0,1)
ax.set_xlim(1,10)
ax.tick_params(axis='x', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(autominorLocator())
ax.tick_params(axis='y', labelsize=32,width=4,colors='black',direction="in",grid_color='black',which='major', length=9)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y',which='minor', length=5,width=4,direction='in')
ax.tick_params(axis='x',which='minor', length=5,width=4,direction='in')
ax.xaxis.grid(False, which='minor')
#ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_ticks(np.arange(1, 11, 1))
#ax.text(0.,0.65,'a)',fontsize=32,fontweight='bold')
legend=plt.legend(loc='best',fontsize=24,handlelength=2, handletextpad=1)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_edgecolor('black')
for text in legend.get_texts():
    text.set_color('black')
    text.set_weight('bold')
plt.show()
