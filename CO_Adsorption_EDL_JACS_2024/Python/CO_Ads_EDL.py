# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 00:33:07 2023

Title:Investigating potential-dependent CO* adsorption using analytical GC-DFT. 
Requirement: The "DFT_CO_data.xlsv" script is used as a template for extracting data through pandas. Feel free to modify the code as is and make sure the paths are corrected.

This script extract data from an excel sheet to create the following plots:
1. Potential-dependent energy change w.r.t one set of EDL properties: dielectric constant and EDL width
2. Potential-independent raw adsorption energies of CO*
3. Compartmentalization of each EDL effect term across different models
4. Sensitivity of free energy change w.r.t potential for a range of EDL properties
5. Main text figure: Potential-dependent CO* adsorption acorss different EDL and adsorption path models

Note: A presummed Helmholtz model is used due to its simplicity. Feel free to rederive in terms of an EDL model of interest.

Feel free to reach out if you have any concerns.
@author: Andrew Jark-Wah Wong (Email: ajwongphd@gmail.com)
"""
# Import Packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#%% Import Data

sheet = 'CoAds.py'
path = "C:/Users/Drew Wong/OneDrive/Documents/GitHub/dft-electrocatalysis-Janik/CO_Adsorption_EDL_JACS_2024/Excel_Data_Reference/DFT_CO_data.xlsx"
#Import Data
df = pd.read_excel(path, sheet_name=sheet)


M= df['Final'].tolist() #Labels

# Initial State Ex: Gas phase CO + Bare Surface (with Li* or Na*)
e_in = df['G_In'].tolist() #Energy of Initial State (eV)
dm_in = df['DM_In'].tolist()  #Dipole moment of Initial state (eA)
polar_in_un = df['Polar_In'].tolist() # Uncorrected Polarizability of Gas phase CO + Bare Surface (with Li* or Na*) (eA^2V-)
polar_bare = polar_in_un[0] # Polarizability of the bare surface (eA^2V-)
polar_in =[i - polar_bare for i in polar_in_un] # Corrected Polarizability of the bare surface (with Li* or Na*) (eA^2V-)



# Final State Ex: Adsorbed State
e_fin = df['G_Fin'].tolist() #Energy of Initial State (eV)
dm_fin = df['DM_Fin'].tolist()  #Dipole moment of Initial state (eA) (eA^2V-)
polar_fin_un = df['Polar_Fin'].tolist() #Uncorrected Polarizability of the CO* (with Li* or Na*)  (eA^2V-)
polar_fin = [i-polar_bare for i in polar_fin_un] # Polarizability of the CO* ( with Li* or Na*)  (eA^2V-)


# Cell Parameters
v = df['Volume'].astype(float)[0] #Volume of Slab in (Angstrom^3)
h = df['Height'].astype(float)[0] #Height of Slab (Angstrom)
wf =df['Work Function'].astype(float)[0] #Work function (eV)
# Cell Parameters
a = v/h # Area of the slab (Angstrom^2)

#Potential Range of Interest
vac_she = 4.6 #Converting from V-Abs to V-SHE
#u_pzc = wf - vac_she
u_pzc =0.462

u_low = -1.5 #Lower Potential (V-SHE)
u_high = 1 #Upper Potential (V_SHE)

u = np.linspace(u_low,u_high,25) # Potential Range (V-SHE)
u_prime = u - u_pzc #Potential range relative to upzc (V-SHE)

# EDL Model Parameters
e_vac = 0.00553 #Vacuum potential
g_solv = 0 #Solvation free energy change of the reaction (eV) 


#Dipole moment and Polarizability Changes

diff_dm = [i - y for i,y in zip(dm_fin,dm_in)]
diff_dm_sq = [i**2 - y**2 for i,y in zip(dm_fin,dm_in)]
diff_polar = [i - y for i,y in zip(polar_fin,polar_in)]
diff_dm_polar_sq = [y * x * x - z * a * a for x,y,z,a in zip(dm_fin,polar_fin,dm_in,polar_in)]
diff_a_dm = [x * y - z* a for x,y,z,a in zip(dm_fin,polar_fin,dm_in,polar_in)]

#%% Solve analytical GC-DFT Math based on one set of er and d

#Calculate aGC-DFT equations
results = {}
# Specify the combinations of er_val and d_val for each M_val
combinations = {
    0: [ (1, 6)],
    1: [ (1, 6)],
    2: [ (1, 6)]
}


# Nested loops to iterate over combinations of er_val, d_val, and M_val
for i, combo_list in combinations.items():
    for er_val, d_val in combo_list:
        # Calculate and store results only for specific combinations of er_val and d_val
        g_1a = e_fin[i] - e_in[i] + g_solv 
        
        # Calculate g_1b for the current iteration
        g_1b = g_1a 
        
        # Calculate capacitance terms
        C = er_val * a * e_vac / d_val
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
        p_0_d = 2 * ((er_val*e_vac) ** 2) * (a ** 2) * (d_val** 2)
        p_0 = diff_dm_polar_sq / p_0_d
        # Polarizability 1st order
        p_1_d = er_val*e_vac * a * d_val ** 2
        p_1 = -1*u_prime * diff_a_dm / p_1_d
        # Polarizability 2nd order
        p_2 = 0.5 * (u_prime ** 2) * (diff_polar[i]) * (1 / d_val) * (1 / d_val)
        # Polarizability total
        p_total = p_0 + p_1 + p_2
        # Calculate G_2C
        g_2c = g_2b + p_total
        EDL_total = p_total + dm_total+c_total

        # Store the result in the dictionary
        results[M[i]] = {'u': u,'g_1a':g_1a,'g_1b': g_1b, 'g_2c': g_2c,'c_total':c_total,'dm_total':dm_total,'p_total':p_total,'EDL_total':EDL_total}


        
#%% Plot Potential-depndent adsorptioin
# Plotting
colors = sns.color_palette('deep', n_colors=30)
fig = plt.figure(figsize=(18, 16))
ap2 = fig.add_subplot(111) #rows,columns,subgraph

# Enumerate over M indices
for i, (M, data) in enumerate(results.items()):
    # Plot adsorption energy w.r.t potential for each g_2c aka full work up EDL value
    if i == 0: #used if and else for colors sequencing based on seaborns
        plt.plot(data['u'], data['g_2c'],
             label=f'{M}',
             linewidth=9,color=colors[i],linestyle='-')      
    else: 
        plt.plot(data['u'], data['g_2c'],
             label=f'{M}',
             linewidth=9,color=colors[i+1],linestyle='-')
# For Clarity H line
plt.hlines(y=0, xmin=-1.5, xmax=0.25, linewidth=5, alpha=0.6, linestyle="--", color='grey')

# Set labels and legend
plt.xlabel("U (V-SHE)", fontsize=48)
plt.ylabel("Adsorption Energy (eV)", fontsize=48)
plt.legend(fontsize=36,ncol=1, loc='best')


plt.text(-1.43, .225, 'a)', fontsize=45,fontweight='bold')

plt.ylim([-.6, 0.3])
plt.xlim([-1.5, 0.2])
ap2.xaxis.set_ticks(np.arange(-1.5, 0.5, 0.25))



# Set axis parameters
ap2.tick_params(axis='x', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=24)
ap2.tick_params(axis='y', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=24)
ap2.tick_params(axis='y', which='minor', length=3, width=2, direction='in')
ap2.tick_params(axis='x', which='minor', length=3, width=2, direction='in')
ap2.set_facecolor('white')
ap2.patch.set_edgecolor('black')
ap2.patch.set_linewidth(5)
# Set edge color of spines to black
for spine in ap2.spines.values():
    spine.set_edgecolor('black')
# Adjust layout and display plot
#plt.tight_layout()
plt.grid(False)
plt.show()

#%% Potential-independent raw adsorption energy values of CO* 

# Plotting
colors = sns.color_palette('deep', n_colors=30)
fig = plt.figure(figsize=(18, 16))
ap2 = fig.add_subplot(111) #rows,columns,subgraph


# Enumerate over M indices
for i, (M, data) in enumerate(results.items()):
    # Plot hlines for each g_2b value
    if M == 'CO*':
        plt.hlines(y=data['g_1b'], xmin=-1.5, xmax=0.25,label=f'{M}', linewidth=8, color=colors[i], linestyle='-')
    else:
        plt.hlines(y=data['g_1b'], xmin=-1.5, xmax=0.25,label=f'{M}', linewidth=8, color=colors[i+1], linestyle='-')
#H lines clarity
plt.hlines(y=0, xmin=-1.5, xmax=0.25, linewidth=7, alpha=1, linestyle='--', color='grey')

# Set labels and legend
plt.xlabel("U (V-SHE)", fontsize=48)
plt.ylabel("Adsorption Energy (eV)", fontsize=48)
plt.legend(fontsize=36,ncol=1, loc='best')


plt.hlines(y=0, xmin=-1.5, xmax=0.25, linewidth=7, alpha=1, linestyle='--', color='grey')
plt.ylim([-.3, 0.2])
plt.xlim([-1.5, 0.2])
ap2.yaxis.set_ticks(np.arange(-0.5, 0.3, 0.1))
ap2.xaxis.set_ticks(np.arange(-1.5, 0.5, 0.25))
# Add text for er and d
#plt.text(-1.4, .85, f'$\epsilon_r$ = {er_values[0]}', fontsize=20)
#plt.text(-1.4, 0.7, f'd = {d_values[0]} $\AA$', fontsize=20)

# Set axis parameters
ap2.tick_params(axis='x', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=25)
ap2.tick_params(axis='y', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=25)
ap2.tick_params(axis='y', which='minor', length=3, width=2, direction='in')
ap2.tick_params(axis='x', which='minor', length=3, width=2, direction='in')
ap2.set_facecolor('white')
ap2.patch.set_edgecolor('black')
ap2.patch.set_linewidth(5)

# Adjust layout and display plot
#plt.tight_layout()
plt.grid(False)
plt.show()


#%% Compartmentalization of EDL effects

# Plotting
colors = sns.color_palette('deep', n_colors=30)
fig, axes = plt.subplots(nrows=len(results), ncols=3, figsize=(18, 14), sharex=True)

#Plots the capacitance, dipole-field, and polarizability contribution for each model w.r.t potential
for i, ((M), data) in enumerate(results.items()):
    # If and else for colors alignment with seaborns but code is the same
    if i==0:
        axes[i, 0].plot(data['u'], data['c_total'], label='c_total', linewidth=7, color=colors[i])
        axes[i, 0].set_ylabel(r"$\Delta$G Capacitive (eV)", fontsize=26)
        axes[i,0].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
        axes[i,0].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)  
        # Plot dm_total
        axes[i, 1].plot(data['u'], data['dm_total'], label='dm_total', linewidth=7, color=colors[i])
        axes[i, 1].set_ylabel(r"$\Delta$G Dipole-Field (eV)", fontsize=26)
        axes[i,1].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
        axes[i,1].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)    
        # Plot p_total
        axes[i, 2].plot(data['u'], data['p_total'], label='p_total', linewidth=7, color=colors[i])
        axes[i, 2].set_ylabel(r"$\Delta$G Polarizability (eV)", fontsize=25)
        axes[i,2].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
        axes[i,2].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
    else:
        axes[i, 0].plot(data['u'], data['c_total'], label='c_total', linewidth=7, color=colors[i+1])
        axes[i, 0].set_ylabel(r"$\Delta$G Capacitive (eV)", fontsize=26)
        axes[i,0].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
        axes[i,0].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)  
        # Plot dm_total
        axes[i, 1].plot(data['u'], data['dm_total'], label='dm_total', linewidth=7, color=colors[i+1])
        axes[i, 1].set_ylabel(r"$\Delta$G Dipole-Field (eV)", fontsize=26)
        axes[i,1].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
        axes[i,1].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)    
        # Plot p_total
        axes[i, 2].plot(data['u'], data['p_total'], label='p_total', linewidth=7, color=colors[i+1])
        axes[i, 2].set_ylabel(r"$\Delta$G Polarizability (eV)", fontsize=25)
        axes[i,2].tick_params(axis='y', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)
        axes[i,2].tick_params(axis='x', labelsize=24,width=4,colors='black',direction="in",grid_color='black',which='major', length=9,pad=10)

# Set common labels and parameters
for ax in axes[-1, :]:
    ax.set_xlabel("U (V-SHE)", fontsize=28)
    ax.set_xlim([-1.5, 0.45])
    ax.xaxis.set_ticks(np.arange(-1.5, 0.25, 0.5))

# Set y-axis limits for all subplots
for ax_row in axes:
    for ax in ax_row:
        ax.set_ylim([-0.1, 0.3])
        ax.set_facecolor('white')
        ax.patch.set_edgecolor('black')
        ax.patch.set_linewidth(5)


# Set layout and display
plt.tight_layout()
plt.show()

#%% Sensitivity of a particular adsorption model due to the EDL properties

M= df['Final'].tolist() #Labels
# Select different dielectric constants and EDL widths
er = [1,2,8,13,78.4] #Relative permittivity (Dielectric Constant)
d = [3,4.5,6,1000] #Helmholtz EDL Width in Angstrom

# Define lists of values for er and d
er_values = er  # Relative permittivity (Dielectric Constant)
d_values = d  # Helmholtz EDL Width in Angstrom

# Initialize results dictionary
results = {}

# Nested loops to iterate over er and d values
for er_val in er_values:
    for d_val in d_values:
        for i in range(len(M)):
            # Calculate g_1a for the current iteration
            g_1a = e_fin[i] - e_in[i] + g_solv 
            
            # Calculate g_1b for the current iteration
            g_1b = g_1a 
            
            # Calculate capacitance terms
            C = er_val * a * e_vac / d_val
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
            p_0_d = 2 * ((er_val*e_vac) ** 2) * (a ** 2) * (d_val** 2)
            p_0 = diff_dm_polar_sq / p_0_d
            # Polarizability 1st order
            p_1_d = er_val*e_vac * a * d_val ** 2
            p_1 = -1*u_prime * diff_a_dm / p_1_d
            # Polarizability 2nd order
            p_2 = 0.5 * (u_prime ** 2) * (diff_polar[i]) * (1 / d_val) * (1 / d_val)
            # Polarizability total
            p_total = p_0 + p_1 + p_2
            # Calculate G_2C
            g_2c = g_2b + p_total
            EDL_total = p_total + dm_total+c_total

            # Store the result in the dictionary
            results[(er_val, d_val, M[i])] = {'u': u, 'g_2c': g_2c, 'c_total': c_total, 'dm_total': dm_total, 'p_total':p_total,'EDL_total': EDL_total}


# Plotting



colors = sns.color_palette('deep', n_colors=30)
fig = plt.figure(figsize=(18, 16))
ap2 = fig.add_subplot(111) #rows,columns,subgraph

import matplotlib.pyplot as plt
import seaborn as sns

# Your existing code

# Define line styles and colors
line_styles = ['-', '--', '-.', ':']  # Different line styles
colors = sns.color_palette('deep', n_colors=30)  # Get a palette of colors

# Initialize counters for line styles and colors
line_style_counter = 0
color_counter = 0

# Filter results for a specific M[i] (let's say M[0])
M_to_plot = M[0]
filtered_results = {key: value for key, value in results.items() if key[2] == M_to_plot}

# Enumerate over combinations of er_val and d_val
for i, ((er_val, d_val, M_val), data) in enumerate(filtered_results.items()):
    # Get the current line style and color
    line_style = line_styles[line_style_counter]
    color = colors[color_counter]

    # Plot u against g_2b with the current line style and color
    plt.plot(data['u'], data['g_2c'],
             label=fr"e$_r$={er_val}, d={d_val} Å",
             linewidth=5,
             linestyle=line_style,
             color=color)

    # Increment counters for line styles and colors
    line_style_counter = (line_style_counter + 1) % len(line_styles)
    color_counter = (color_counter + 1) % len(colors)

# Set labels and legend
plt.xlabel("U (V-SHE)", fontweight='bold', fontsize=48)
plt.ylabel("Adsorption Energy (eV)", fontweight='bold', fontsize=48)
plt.legend(fontsize=25,ncol=3, loc='best')


plt.hlines(y=0, xmin=-1.5, xmax=0.25, linewidth=7, alpha=1, linestyle='--', color='grey')
plt.ylim([-.6, 0.4])
plt.xlim([-1.5, 0.2])
ap2.yaxis.set_ticks(np.arange(-0.6, 0.6, 0.1))
ap2.xaxis.set_ticks(np.arange(-1.5, 0.5, 0.25))
# Add text for er and d
#plt.text(-1.4, .85, f'$\epsilon_r$ = {er_values[0]}', fontsize=20)
#plt.text(-1.4, 0.7, f'd = {d_values[0]} $\AA$', fontsize=20)

# Set axis parameters
ap2.tick_params(axis='x', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=25)
ap2.tick_params(axis='y', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=25)
ap2.tick_params(axis='y', which='minor', length=3, width=2, direction='in')
ap2.tick_params(axis='x', which='minor', length=3, width=2, direction='in')
ap2.set_facecolor('white')
ap2.patch.set_edgecolor('black')
ap2.patch.set_linewidth(8)

# Adjust layout and display plot
#plt.tight_layout()
plt.grid(False)
plt.show()

# %% Main text figure: Potential-dependent CO* adsorption across different models of EDL and adsorption path

#Initialize empty dictionary
results = {}

colors = sns.color_palette('deep', n_colors=30)
# Specify the combinations of er_val and d_val for each M_val
combinations = {
    0: [(78.4, 3), (1, 6)],
    1: [(78.4, 3), (1, 6)],
    2: [(78.4, 3), (1, 6)]
}

# Nested loops to iterate over combinations of er_val, d_val, and M_val
for i, combo_list in combinations.items():
    for er_val, d_val in combo_list:
        # Calculate and store results only for specific combinations of er_val and d_val
        g_1a = e_fin[i] - e_in[i] + g_solv 
        
        # Calculate g_1b for the current iteration
        g_1b = g_1a 
        
        # Calculate capacitance terms
        C = er_val * a * e_vac / d_val
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
        p_0_d = 2 * ((er_val*e_vac) ** 2) * (a ** 2) * (d_val** 2)
        p_0 = diff_dm_polar_sq / p_0_d
        # Polarizability 1st order
        p_1_d = er_val*e_vac * a * d_val ** 2
        p_1 = -1*u_prime * diff_a_dm / p_1_d
        # Polarizability 2nd order
        p_2 = 0.5 * (u_prime ** 2) * (diff_polar[i]) * (1 / d_val) * (1 / d_val)
        # Polarizability total
        p_total = p_0 + p_1 + p_2
        # Calculate G_2C
        g_2c = g_2b + p_total
        EDL_total = p_total + dm_total+c_total

        # Store the result in the dictionary
        results[(er_val, d_val, M[i])] = {'u': u, 'g_2c': g_2c, 'c_total': c_total, 'dm_total': dm_total, 'p_total':p_total,'EDL_total': EDL_total}



# Define colors for each M_val

# Plotting
fig = plt.figure(figsize=(22, 16))
ax = fig.add_subplot(111)
M_val = range(len(M))

lw = 7
# Keep track of which M_names have been plotted
plotted_M_names = set()

# Iterate through results and plot each line
for key, data in results.items():
    er_val, d_val, M_name = key
    if M_name not in plotted_M_names:
        # Plot the line and add it to the plotted_M_names set
        if M_name == M[0]:
            line1, = ax.plot(data['u'], data['g_2c'], color=colors[0], linewidth=lw, linestyle='dashed', label=M_name + fr": $\epsilon$$_r$={er_val}, d={d_val} Å",alpha=0.6)
        elif M_name == M[1]:
            line2, = ax.plot(data['u'], data['g_2c'], color=colors[2], linewidth=lw, linestyle='dashed', label=M_name + fr": $\epsilon$$_r$={er_val}, d={d_val} Å",alpha=0.6)
        else:
            line3, = ax.plot(data['u'], data['g_2c'], color=colors[3], linewidth=lw, linestyle='dashed', label=M_name + fr": $\epsilon$$_r$={er_val}, d={d_val} Å",alpha=0.6)

        plotted_M_names.add(M_name)
    else:
        # Plot the line without adding to the legend
        if M_name == M[0]:
            ax.plot(data['u'], data['g_2c'], linewidth=lw, linestyle='solid', color=colors[0],label=M_name + fr": $\epsilon$$_r$={er_val}, d={d_val} Å")
            #ax.fill_between(data['u'], line1.get_ydata(), data['g_2c'], color=colors[0], alpha=0.2)
        elif M_name == M[1]:
            ax.plot(data['u'], data['g_2c'], linewidth=lw, linestyle='solid', color=colors[2],label=M_name + fr": $\epsilon$$_r$={er_val}, d={d_val} Å")
            #ax.fill_between(data['u'], line2.get_ydata(), data['g_2c'], color=colors[1], alpha=0.2)
        else:
            ax.plot(data['u'], data['g_2c'], linewidth=lw, linestyle='solid', color=colors[3],label=M_name + fr": $\epsilon$$_r$={er_val}, d={d_val} Å")
            #ax.fill_between(data['u'], line3.get_ydata(), data['g_2c'], color=colors[2], alpha=0.2)



# Set labels and legend
plt.xlabel("U (V-SHE)", fontsize=48)
plt.ylabel("Adsorption Energy (eV)", fontsize=48)
plt.legend(fontsize=28,ncol=2, loc='best',fancybox=True,facecolor='white')


plt.hlines(y=0, xmin=-1.5, xmax=0.25, linewidth=5, alpha=0.6, linestyle="--", color='grey')
plt.ylim([-.6, 0.3])
plt.xlim([-1.5, 0.2])
ax.yaxis.set_ticks(np.arange(-0.6, 0.4, 0.1))
ax.xaxis.set_ticks(np.arange(-1.5, 0.5, 0.25))


# Set axis parameters
ax.tick_params(axis='x', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=25)
ax.tick_params(axis='y', labelsize=44, width=5, colors='black', direction="in", grid_color='black', which='major', length=10, pad=25)
ax.tick_params(axis='y', which='minor', length=3, width=2, direction='in')
ax.tick_params(axis='x', which='minor', length=3, width=2, direction='in')
ax.set_facecolor('white')
ax.patch.set_edgecolor('black')
ax.patch.set_linewidth(5)

#plt.tight_layout()
plt.grid(False)
plt.show()
