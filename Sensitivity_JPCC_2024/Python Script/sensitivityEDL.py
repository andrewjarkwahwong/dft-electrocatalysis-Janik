# -*- coding: utf-8 -*-
"""
Created on Friday Feb 15 2024

Analytical GC framework for CO reduction on Cu 111/Cu 100

Using the Excel Template, this code will provide 3 figures based on the analytical GC-DFT approach:
1. Compartmentalization of potential-dependent contribution on the barrier
2. Sensitivity of barriers based on EDL properties
3. Sensitivity of symmetry factors based on EDL Properties

Usage of script and the analytical GC-DFT framework requires citations of both work:
1. https://doi.org/10.1016/j.jcat.2024.115360
2. (JPCC Once it's published)

@author: Andrew Jark-Wah Wong (Email: ajwongphd@gmail.com)
"""
# %% User Inputs
import PySimpleGUI as sg
import ipywidgets 
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate as tb
import pandas as pd
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator

# Sheet name needs to be specified
sheet = '111.py' # 100.py or 111.py
#Import Data
df = pd.read_excel("C:/Users/Drew Wong/OneDrive/Documents/GitHub/dft-electrocatalysis-Janik/Sensitivity_JPCC_2024/Excel Sheet/CO_data.xlsx", sheet_name=sheet)


M = df['M'].tolist() #Reaction Names

# Bare Metals
# Area and Upzc 
a=  df['Area'].iloc[0] #50.79 for 111 and 58.64 for 100
u_pzc = df['Upzc'].iloc[0] #0.29 and -0.15 for 111 and 100 
polar_bare = df['Polar_Bare'].iloc[0] #Polarizability of the Bare Surface

# Initial State Ex: CO* + H2 (g) for OC-H formation
e_in = df['E_In'].tolist() #Energy of Initial State (eV)
dm_in = df['DM_In'].tolist()  #Dipole moment of Initial state (eA)
polar_in_un = df['Polar_In'].tolist() #Polarizability of Initial state (eA^2V^-1)
polar_in =[i - polar_bare for i in polar_in_un] #Polarizability of the adsorbates

# Transition State
e_fin = df['E_Fin'].tolist() #Energy of Transition State (eV)
dm_fin = df['DM_Fin'].tolist()  #Dipole moment of Transition State (eA)
polar_fin_un = df['Polar_Fin'].tolist() #Polarizability of the Transition State
polar_fin = [i-polar_bare for i in polar_fin_un] #Polarizability of the adsorbates

g_solv = df['G_Solv'].tolist() #solvation free energy change of the reaction (eV)



#Potential Range of Interest and general inputs
vac_nhe = 4.6 #Converting from V-Abs to V-NHE
u_low = -2.5 #Lower Potential V-NHE
u_high = 1 #Upper Potential V-NHE

# EDL Model Parameters
e_vac = 0.00553 #Vacuum permittivity 


# %% Figure 3: Compartmentalization (The Math)

#USER INPUT: What dielectric constants and d for for each reaction
M_values = {'C-H': {'er': 78.4, 'd': 3}, 'O-H': {'er': 78.4, 'd': 3}, 'OC-CO': {'er': 78.4, 'd': 3}}

# Run the rest of the code and generate the math
#Dipole moment and Polarizability Changes
diff_dm = [i - y for i,y in zip(dm_fin,dm_in)]
diff_polar = [i - y for i,y in zip(polar_fin,polar_in)]
diff_dm_polar_sq = [y * x * x - z * a * a for x,y,z,a in zip(dm_fin,polar_fin,dm_in,polar_in)]
diff_dm_sq = [i ** 2 - y ** 2 for i,y in zip(dm_fin,dm_in)] 
diff_a_dm = [x * y - z* a for x,y,z,a in zip(dm_fin,polar_fin,dm_in,polar_in)]
             
#Potential 
u = np.linspace(u_low,u_high,150)
u_prime = u - u_pzc 

#Store into a dictionary called results
results = {}

# Solve for G_2c vs U based on different er and d combinations (Should make this a function)
for i in range(len(M)):
    # Get parameters for er and d based on M
    er = M_values[M[i]]['er']
    d = M_values[M[i]]['d']

    if i == 2:  # For C-C chemical Step (non-faradaic)

        # Model 1A: Free Energy Change at U_PZC 
        g_1a = e_fin[i] - e_in[i] + g_solv[i] 
        
        # Model 1B: Potenital dependent Free Energy Change assuming Beta = 0 (no e- transfer)
        g_1b = g_1a 
        
        #Model 2A: Incorporate Capacitance Charge
        C = er*a*e_vac/d # Predicted Capacitiance by Helmholtz Model
        diff_dm_sq = (dm_fin[i] ** 2 - dm_in[i] ** 2) 
        # 0th order capacitance
        C_0 = -0.5 * (diff_dm_sq) / (C * d ** 2)
        #  1st order wrt U'
        C_const_1 = diff_dm[i] / d
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
        diff_dm_polar_sq = (polar_fin[i] * dm_fin[i] * dm_fin[i] - polar_in[i] * dm_in[i] * dm_in[i])
        diff_a_dm = polar_fin[i] * dm_fin[i] - polar_in[i] * dm_in[i]
        # Polarizability 0th order
        p_0_d = 2 * ((er*e_vac) ** 2) * (a ** 2) * (d** 2)
        p_0 = diff_dm_polar_sq / p_0_d
        # Polarizability 1st order
        p_1_d = er*e_vac * a * d ** 2
        p_1 = -1*u_prime * diff_a_dm / p_1_d
        # Polarizability 2nd order
        p_2 = 0.5 * (u_prime ** 2) * (diff_polar[i]) * (1 / d) * (1 / d)
        # Polarizability total
        p_total = p_0 + p_1 + p_2
        # Calculate G_2C
        g_2c = g_2b + p_total
        EDL_total = p_total + dm_total+c_total
        
        # Now you can use g_2c for further analysis or store the results as needed
        # Store the result in the dictionary
        results[M[i],er,d] = {'u': u, 'g_1a':g_1a,'g_1b': g_1b, 'g_2c': g_2c,
                          'c_total': c_total, 'dm_total': dm_total, 'p_total': p_total, 'EDL_total': EDL_total}
        
    else:  # For C-H and O-H Step (Faradaic Steps)
        # Model 1A: Free Energy Change at U_PZC
        g_1a = e_fin[i] - e_in[i] + g_solv[i] + u_pzc
        
        # Model 1B: Potenital dependent Free Energy Change assuming Beta = 1 (Full e- transfer)
        g_1b = g_1a+u_prime
        
        #Model 2A: Incorporate Capacitance Charge
        C = er*a*e_vac/d # Predicted Capacitiance by Helmholtz Model
        diff_dm_sq = (dm_fin[i] ** 2 - dm_in[i] ** 2) 
        # 0th order capacitance
        C_0 = -0.5 * (diff_dm_sq) / (C * d ** 2)
        #  1st order wrt U'
        C_const_1 = diff_dm[i] / d
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
        diff_dm_polar_sq = (polar_fin[i] * dm_fin[i] * dm_fin[i] - polar_in[i] * dm_in[i] * dm_in[i])
        diff_a_dm = polar_fin[i] * dm_fin[i] - polar_in[i] * dm_in[i]
        # Polarizability 0th order
        p_0_d = 2 * ((er*e_vac) ** 2) * (a ** 2) * (d** 2)
        p_0 = diff_dm_polar_sq / p_0_d
        # Polarizability 1st order
        p_1_d = er*e_vac * a * d ** 2
        p_1 = -1*u_prime * diff_a_dm / p_1_d
        # Polarizability 2nd order
        p_2 = 0.5 * (u_prime ** 2) * (diff_polar[i]) * (1 / d) * (1 / d)
        # Polarizability total
        p_total = p_0 + p_1 + p_2
        # Calculate G_2C
        g_2c = g_2b + p_total
        EDL_total = p_total + dm_total+c_total
        
        # Now you can use g_2c for further analysis or store the results as needed
        # Store the result in the dictionary along with er and d
        results[M[i],er,d] = {'u': u, 'g_1a':g_1a,'g_1b': g_1b, 'g_2c': g_2c,
                          'c_total': c_total, 'dm_total': dm_total, 'p_total': p_total, 'EDL_total': EDL_total}

# %% Figure 3: Compartmentalization (The Plot)

#USER INPUT: Change Voltage Here
volts=[-0.5,-0.75] 


#Run the cell and generate the plot
u = np.linspace(u_low,u_high,150)
indexv = [np.argmin(np.abs(u - volt)) for volt in volts]

u_prime_1=[volts[0]-u_pzc,volts[0]-u_pzc,0]
u_prime_2=[volts[1]-u_pzc,volts[1]-u_pzc,0]

# Extract values from 
g_1a_values_1 = [results[(M, er, d)]['g_1a'] for M in M_values]
c_total_values_1 = [results[(M, er, d)]['c_total'][indexv][0] for M in M_values]
dm_total_values_1 = [results[(M, er, d)]['dm_total'][indexv][0] for M in M_values]
p_total_values_1 = [results[(M, er, d)]['p_total'][indexv][0] for M in M_values]
g_2c_values_1 = [results[(M, er, d)]['g_2c'][indexv][0] for M in M_values]

g_1a_values_2 = [results[(M, er, d)]['g_1a'] for M in M_values]
c_total_values_2 = [results[(M, er, d)]['c_total'][indexv][1] for M in M_values]
dm_total_values_2 = [results[(M, er, d)]['dm_total'][indexv][1] for M in M_values]
p_total_values_2 = [results[(M, er, d)]['p_total'][indexv][1] for M in M_values]
g_2c_values_2 = [results[(M, er, d)]['g_2c'][indexv][1] for M in M_values]

# Plotting
bar_width = 0.35
index = np.arange(len(M_values))
colors = sns.color_palette('deep', n_colors=9)
fig, ax = plt.subplots(figsize=(24, 19))

# Plot for -0.5 V-SHE
bar1 = ax.bar(index+ bar_width, g_1a_values_1, bar_width,  color=colors[1], alpha=0.5)
bar5 = ax.bar(index+ bar_width, u_prime_1, bar_width,  color=colors[2], alpha=0.5)
ax.axhline(y=0, xmin=0, xmax=10, linewidth=5, alpha=0.9, color='black', linestyle='--')
ax.scatter(index+ bar_width, g_2c_values_1, s=1000, label='Total Barrier at'+ str(volts[0]) +'V-SHE', color='white', marker='o',edgecolor='black',linewidth=7)

# Plot for -0.75 V-SHE
bar1_7 = ax.bar(index, g_1a_values_2, bar_width, label=r'$\Delta$G$^{‡o}$',color=colors[1], alpha=1)
bar5_7 = ax.bar(index, u_prime_2, bar_width, color=colors[2],label='|e|U$^{\prime}$', alpha=1)
ax.scatter(index, g_2c_values_2,label='Total Barrier at'+ str(volts[1]) +'V-SHE', s=1000, color='white', marker='D', alpha=1,edgecolor='black',linewidth=7,zorder=3200)

# Initialize legend_added flag
legend_added = False

for i in range(len(M_values)):
    if i == 1:
        # Stack bars differently for index O-H
        bar2 = ax.bar(index[i]+ bar_width, c_total_values_1[i], bar_width,  bottom=g_1a_values_1[i], color=colors[3], alpha=0.5)
        bar3 = ax.bar(index[i]+ bar_width, dm_total_values_1[i], bar_width,  bottom=g_1a_values_1[i] + c_total_values_1[i], color=colors[0], alpha=0.5)
        bar4 = ax.bar(index[i]+ bar_width, p_total_values_1[i], bar_width,  bottom=g_1a_values_1[i] + c_total_values_1[i] + dm_total_values_1[i], color=colors[6], alpha=0.5)
        ax.text(index[i]-0.5+ bar_width,g_1a_values_1[i] + c_total_values_1[i] + dm_total_values_1[i]+0.25,'-0.5 V',fontsize=36,fontweight='bold')
        ax.text(index[i]+0.225,g_1a_values_2[i] + c_total_values_2[i] + dm_total_values_2[i]-0.1,'-0.25 V',fontsize=36,fontweight='bold')
        # Bars for Index 7 (shifted by bar width)
        bar2_2 = ax.bar(index[i], c_total_values_2[i], bar_width, bottom=g_1a_values_2[i], color=colors[3],label='Capacitive', alpha=1)
        bar3_2 = ax.bar(index[i], dm_total_values_2[i], bar_width,  bottom=g_1a_values_2[i] + c_total_values_2[i], color=colors[0],label='Dipole-Field', alpha=1)
        bar4_2 = ax.bar(index[i], p_total_values_2[i], bar_width,  bottom=g_1a_values_2[i] + c_total_values_2[i] + dm_total_values_2[i], label='Polarizability',color=colors[6], alpha=1)
    elif i == 0:
        # Stack bars differently for other indices for C-H
        bar2 = ax.bar(index[i]+ bar_width, c_total_values_1[i], bar_width, bottom=u_prime_1[i], color=colors[3], alpha=0.5)
        bar3 = ax.bar(index[i]+ bar_width, dm_total_values_1[i], bar_width, bottom=u_prime_1[i] + c_total_values_1[i], color=colors[0], alpha=0.5)
        bar4 = ax.bar(index[i]+ bar_width, p_total_values_1[i], bar_width, bottom=u_prime_1[i] + c_total_values_1[i] + dm_total_values_1[i], color=colors[6], alpha=0.5)
        ax.text(index[i]-0.5+ bar_width,g_1a_values_1[0]+0.1,'-0.5 V',fontsize=36,fontweight='bold')
        ax.text(index[i]+0.225,g_1a_values_1[0]+0.1,'-0.25 V',fontsize=36,fontweight='bold')
        # Bars for Index 7 (shifted by bar width)
        bar2_2 = ax.bar(index[i], c_total_values_2[i], bar_width, bottom=u_prime_2[i], color=colors[3], alpha=1)
        bar3_2 = ax.bar(index[i], dm_total_values_2[i], bar_width, bottom=u_prime_2[i] + c_total_values_2[i], color=colors[0], alpha=1)
        bar4_2 = ax.bar(index[i], p_total_values_2[i], bar_width, bottom=u_prime_2[i] + c_total_values_2[i] + dm_total_values_2[i], color=colors[6], alpha=1)
    else: #C-C
        bar2 = ax.bar(index[i]+ bar_width, c_total_values_1[i], bar_width, color=colors[3], alpha=0.5)
        bar3 = ax.bar(index[i]+ bar_width, dm_total_values_1[i], bar_width, bottom=c_total_values_1[i], color=colors[0], alpha=0.5)
        bar4 = ax.bar(index[i]+ bar_width, p_total_values_1[i], bar_width, bottom=c_total_values_1[i] + dm_total_values_1[i], color=colors[6], alpha=0.5)
        ax.text(index[i]-0.5+ bar_width,g_1a_values_1[2]+0.1,'-0.5 V',fontsize=36,fontweight='bold')
        ax.text(index[i]+0.225,g_1a_values_1[2]+0.1,'-0.25 V',fontsize=36,fontweight='bold')
        # Bars for Index 7 (shifted by bar width)
        bar2_2 = ax.bar(index[i], c_total_values_2[i], bar_width, color=colors[3], alpha=1)
        bar3_2 = ax.bar(index[i], dm_total_values_2[i], bar_width, bottom=c_total_values_2[i], color=colors[0], alpha=1)
        bar4_2 = ax.bar(index[i], p_total_values_2[i], bar_width, bottom=c_total_values_2[i] + dm_total_values_2[i], color=colors[6], alpha=1)

    # Add legend only if not added before
    if not legend_added:
        ax.legend()
        legend_added = True

# Customize the plot
#ax.set_xlabel('M[i]')
ax.set_ylabel('Free Energy (eV)',fontsize=52, fontweight ='bold')
#ax.set_title('Stacked Bar Chart of Energy Components')
ax.set_xticks(index)
ax.set_xticklabels(M_values, fontweight ='bold')

ax.yaxis.set_ticks(np.arange(-2, 2.5, 0.5))
ax.set_ylim([-2, 2.25])
# Set minor ticks on the y-axis
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

ax.set_facecolor('white')
ax.patch.set_edgecolor('black')
ax.patch.set_linewidth(5)
ax.tick_params(axis='both', labelsize=44, width=8, colors='black', direction="in", grid_color='black', which='major', length=20, pad=15)
ax.tick_params(axis='both', which='minor', length=12, width=5, direction='in')
# Set axis parameters for the subplot

# Set ticks on the top and right sides
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

# Set spines visible on all sides
ax.spines['top'].set_linewidth(4)
ax.spines['right'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.spines['left'].set_linewidth(4)


ax.set_facecolor('white')
ax.patch.set_edgecolor('black')
ax.patch.set_linewidth(8)

# Add legend to each subplot with customization
legend= ax.legend(loc='best', fontsize=32, bbox_to_anchor=(1.06,-0.1), borderaxespad=0, ncol=4, frameon=True, fancybox=True, framealpha=1, edgecolor='black', facecolor='white', labelspacing=0.5)
# Customize the legend font weight and box properties
for text in legend.get_texts():
    text.set_fontweight('bold')

legend.get_frame().set_linewidth(3)  # Adjust the legend box linewidth

plt.show()


# %% Figure 4: Sensitivity based on EDL Properties (The Math)

er = [1,2,4,8,13,78.4] #Relative permittivity (Dielectric Constant)
d = [3,4.5,6,10] #Helmholtz EDL Width in Angstrom

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
resultsB = {}
# Nested loops to iterate over er and d values
for er_val in er_values:
    for d_val in d_values:
        for i in range(len(M)):
            if i == 2: #C-C
                # Calculate g_1a for the current iteration
                g_1a = e_fin[i] - e_in[i] + g_solv[i] 
                
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
                
                #beta

                # Beta 
                a_beta = diff_polar[i]/ (2 * d_val ** 2)
                b_beta =  2 * diff_dm[i] / d_val - (diff_a_dm) / (er_val *e_vac* a * d_val ** 2)
                beta = b_beta + 2 * a_beta * u_prime  # eq 31 calculate as function of U wrt to Upzc
                beta_avg = sum(beta) / len(beta)
    
                # Store the result in the dictionary
                resultsB[(er_val, d_val, M[i])] = {'u': u,'g_2c':g_2c, 'beta':beta_avg}
            else: # C-H and O-H
                # Calculate g_1a for the current iteration
                g_1a = e_fin[i] - e_in[i] + g_solv[i] + u_pzc
                
                # Calculate g_1b for the current iteration
                g_1b = g_1a + u_prime
                
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
                
                
                # Beta 
                a_beta = diff_polar[i]/ (2 * d_val ** 2)
                b_beta = 1 + 2 * diff_dm[i] / d_val - (diff_a_dm) / (er_val *e_vac* a * d_val ** 2)
                beta = b_beta + 2 * a_beta * u_prime  # eq 31 calculate as function of U wrt to Upzc
                beta_avg = sum(beta) / len(beta)
                # Store the result in the dictionary
                resultsB[(er_val, d_val, M[i])] = {'u': u,'g_2c':g_2c, 'beta':beta_avg}
# %% Figure 4: Sensitivity based on EDL Properties (The Plot)
# Just run the cell, you may need to adjust d val based on linestyle 
palette = ['Blues','flare','crest']
fig, axs = plt.subplots(1, 3, figsize=(28, 9))
labels=['C-H','O-H','OC-CO']
colors1 = sns.color_palette('deep', n_colors=30)
# Create an empty list to collect the labels for the legend
legend_labels = []
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

for idx, (M_val,diff_dm_val) in enumerate(zip(M[:3],diff_dm[:3])):  # Assuming M is a list
    # Filter results for a specific M[i]
    filtered_results = {key: value for key, value in resultsB.items() if key[2] == M_val}

    # Use a different color palette for each subplot
    colors = sns.color_palette(palette[idx], n_colors=len(filtered_results))

    # Enumerate over combinations of er_val and d_val
    for i, ((er_val, d_val, _), data) in enumerate(filtered_results.items()):
        if d_val == 3:
            linestyle = '-'
        elif d_val == 4.5:
            linestyle = 'dashed'
        elif d_val == 6:
            linestyle = '-.'
        else:
            linestyle = ':'

        # Plot u against g_2b in the current subplot
        line, = axs[idx].plot(data['u'], data['g_2c'],
                              label=fr"e$_r$={er_val}, d={d_val}",
                              linewidth=5, alpha=1, color=colors[i], linestyle=linestyle)
        
        # Append the label to the legend_labels list
        legend_labels.append(fr"e$_r$={er_val}, d={d_val}")

    # Add text in the top right corner of each subplot
    #axs[idx].text(0.08, 0.90, chr(97 + idx)+')', fontsize=30, ha='center', va='center',
                  #transform=axs[idx].transAxes, fontweight='bold')
    if idx == 2:
        axs[idx].text(0.24, 0.92, labels[idx], transform=axs[idx].transAxes, fontsize=36, fontweight='bold', va='center', ha='center', c=colors1[idx])
    else:
        axs[idx].text(0.2, 0.92, labels[idx], transform=axs[idx].transAxes, fontsize=36, fontweight='bold', va='center', ha='center', c=colors1[idx])
    axs[idx].text(0.07, 0.92, f'{chr(97+idx)})', transform=axs[idx].transAxes, fontsize=36, fontweight='bold', va='center', ha='center')

    axs[idx].text(.28, 0.84, fr'$\Delta\mu$ = {diff_dm_val:.2f} '+r'e$\AA$', fontsize=36, ha='center', va='center',
                  transform=axs[idx].transAxes)

    # Set labels for the subplot
    axs[idx].set_xlabel("U (V-SHE)", fontweight='bold', fontsize=32)
    axs[idx].set_ylabel("Activation Free Energy (eV)", fontweight='bold', fontsize=30)

    # Set limits and ticks for the subplot
    axs[idx].set_xlim([-1.5, 0])
    axs[idx].xaxis.set_ticks(np.arange(-1.5, 0.25, 0.25))
    # axs[idx].hlines(y=0, xmin=-1.5, xmax=0, linewidth=2, alpha=0.5, linestyle='--', color='grey')
    if idx == 2: # C-C
        axs[idx].set_ylim([.9, 1.05])
    else:
        axs[idx].set_ylim([0, 1.2])
        #axs[idx].set_ylim([0.2, .6])
    # Set axis parameters for the subplot
    axs[idx].tick_params(axis='both', labelsize=28, width=4, colors='black', direction="in", grid_color='black', which='major', length=10, pad=15)
    axs[idx].tick_params(axis='both', which='minor', length=6, width=5, direction='in')
    
    # Set ticks on the top and right sides
    axs[idx].xaxis.set_ticks_position('both')
    axs[idx].yaxis.set_ticks_position('both')
    
    # Set spines visible on all sides
    axs[idx].spines['top'].set_linewidth(4)
    axs[idx].spines['right'].set_linewidth(4)
    axs[idx].spines['bottom'].set_linewidth(4)
    axs[idx].spines['left'].set_linewidth(4)
    
    axs[idx].tick_params(axis='y', which='minor', length=6, width=5, direction='in')
    axs[idx].tick_params(axis='x', which='minor', length=6, width=5, direction='in')
    
    
    axs[idx].set_facecolor('white')
    minor_locator = AutoMinorLocator(2)
    axs[idx].xaxis.set_minor_locator(minor_locator)
    minor_locator1 = AutoMinorLocator(1)
    axs[idx].yaxis.set_minor_locator(minor_locator1)
    axs[idx].patch.set_edgecolor('black')
    axs[idx].patch.set_linewidth(5)
    # Use a different color palette for each subplot
    colors = sns.color_palette(palette[idx], n_colors=len(filtered_results))


# Adjust layout and display plot
plt.tight_layout()
plt.grid(False)
plt.show()

# %% Symmetry Factor based on EDL Properties
# Specify the desired M[i] based on M (aka the index of the reaction)
desired_M_index = 1  # Change this to the desired index of M
custom_order = [78.4, 13, 8, 4, 2.0, 1.0]  # Replace with your desired order
# Assuming results is a dictionary
data = []
index_values = set()
columns_values = set()

for (er_val, d_val, model), values in resultsB.items():
    if model == M[desired_M_index]:
        index_values.add(er_val)
        columns_values.add(d_val)
        data.append([er_val, d_val, values['beta']])

# Create a DataFrame
ef = pd.DataFrame(data, columns=['er_val', 'd_val', 'beta'])

# Sort the DataFrame based on custom_order
ef['er_val'] = pd.Categorical(ef['er_val'], categories=custom_order, ordered=True)
ef = ef.sort_values(by=['er_val', 'd_val'])

# Pivot the DataFrame for seaborn heatmap
heatmap_data = ef.pivot(index='er_val', columns='d_val', values='beta')


# Change the labels of the index and columns
heatmap_data = heatmap_data.rename_axis('', axis=1)
heatmap_data = heatmap_data.rename_axis('', axis=0)
# Ensure that the annotation values are not rounded too aggressively
heatmap_data_display = heatmap_data.copy()
heatmap_data_display = heatmap_data_display.applymap(lambda x: f"{x:.2f}")
# Create a heatmap with bold labels and changed font size
plt.figure(figsize=(10, 10))
ax = sns.heatmap(heatmap_data, annot=heatmap_data_display, cmap="Spectral", fmt="", linewidths=5,
                 annot_kws={"weight": "bold", "size": 28, "color": "white"}, center=0.5)  # Adjust size as needed

cax = ax.figure.axes[-1]
cax.tick_params(labelsize=24)
ax.tick_params(axis='y', labelsize=24, width=3, colors='black', direction="in", grid_color='black',
               which='major', length=4)
ax.tick_params(axis='x', labelsize=24, width=3, colors='black', direction="in", grid_color='black',
               which='major', length=4)

plt.xlabel('EDL Width (Å)', fontweight='bold', fontsize=24,labelpad=10)
plt.ylabel('Dielectric Constant', fontweight='bold', fontsize=24,labelpad=10)

# Add label to the color bar

cbar = ax.collections[0].colorbar
#cbar.set_label('Average Symmetry Factor', fontsize=24, fontweight='bold', labelpad=20)

# Set color bar limits
ax.collections[0].set_clim(0, 1)

plt.show()