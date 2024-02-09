# Analytical Grand-Canonical DFT Approach for the Calculation of Potential-Dependent Electrochemical Activation Energies

# Background
This part of the repository is to provide a tool to calculate potential dependent electrochemical activation barriers with electrochemical double layer (EDL) consideration using our analytical GC-DFT approach. I have provided several different tools within the repository and links based on your preference, where all tools compute similar sensitivity analysis.

Please read: Our paper in Journal of Catalysis for more details on the theory and derivation of our approach. 
Usage of our approach requires citation of this work. 

The main advantage of our approach is to quantify the sensitivity of DFT prediction activation barriers to different types of EDL models and their parameters. Note that our approach uses a simple Helmholtz model to address both the changes in workfunction along the reaction path and the description of the field (Eqs. 14, 20, and 21). In practice, any model of the EDL, capacitance, and the field can be used and rederived. With a Helmholtz model, we can quantify how different reaction energetics and activation barriers change with the dielectric constant and the EDL width w.r.t potential.


# Available Tools 
Excel Notebook: Excel_Barrier_EDL.xlsx
Jupyter Notebook: Jupyter_Barriers_EDL.ipynb
Python Scripts: Barrier_EDL_Base.py and Barrier_E_and_d.py
Bash scripts: polar and getenergies

--Inputs--

Initial/Final States:
1. DFT Free Energies: Output of DFT energies with corrections (ex: ZPVE, TS) to Free Energies. Note that the reference state and initial state are important to choose. Generally, we use bare surface (*) + $\frac{1}{2}$ H$_$2 and reference to this potential or the H* on the surface and reference to the equilibrium potential of this state. 
2. Dipole Moment: Dipole Moments calculated from VASP (IDIPOL =3, LDIPOL = .TRUE.) for initial and final states. Both should be turned on. I also would test whether
3. Polarizability: Refer to Supplemental Section 3 of how to calculate the polarizability of the surface and the adsorbate

Cell Parameters:
4. Volume: Volume of the bare surface slab (ase gui can provide this quickly via quick info)
5. Height: Total height of the surface slab
6. Fermi Energy: Fermi Energy of the bare surface slab
7. Vacuum Potential: Vacuum Potential determined from either the WAVECAR or CHGCAR (LVHAR=.TRUE.)
8. Bare Metal Polarizability: Polarizability of the bare metal slab with no adsorbate (Same procedure as Input #3)

Potential Range:
9. Lower Limit: Most Negative (reducing) voltage/potential on a NHE scale
10. Upper Limit: Most Positive (oxidative) voltage/potenital on a NHE scale

EDL Model Parameters:
11. Relative Permittivity: Dielectric Constant of the media within the interface ($\epsilon$ $_r$ = 1 for vacuum or $\epsilon$$_r$ = 78.4 for bulk water)
12. Width of the EDL: Width is defined as the distance between the electrode surface and the countercharge ions in angstrom 
13. Vacuum to NHE: The voltage correction from absolute scale to NHE scale (commonly as 4.2 V to 4.8 V but test for your system)
14. Solvation Free Energy Change (eV): This term incorporates the total solvation free energy change along the reaction path. Our approach does not include solvation (micro-solvation was included in our work to model H+ shuttling)


## Excel Notebook: Excel_Barrier_EDL.xlsx
The script provided reproduces the main plots in the manuscript for calculating the activation barrier of NH* to NH$_2$* with 2 H$_2$O molecules. 

--Current Capabilities--
1. Quantification of U$_{pzc}$ of each state along the reaction path using the calculated capacitance of the Helmholtz model
2. Compartmentalized both the potential-dependent and independent EDL term in Figure 3 w.r.t selected $\epsilon_r$ and d. These are quantified for each EDL correction w.r.t potential. 
3. Calculates the sensitivity of the activation barrier w.r.t to potential for different presumed values of beta using model 1b (Figure 4)
4. Compute the finite cell and explicit electrification terms given the $\epsilon_r$ and d (Figure 5)
5. Calculates the barrier profile w.r.t potential between model 1a,2a,2b, and 2c. The slopes are the symmetry factor along the path. (Figure 6)
6. Sensitivity of activation barriers w.r.t potential given different values of $\epsilon_r$ and d (FIgure 10)


--Notes--
1. Shaded Blue are the inputs in the excel sheet. 
2. Work function is calculated as eU$_{vacuum}$ - E_${fermi}$
3. Barriers are calculated using the polarizability w.r.t bare metal, using the polarizability change of only the adsorbate along the rxn path.
4. PZC of each state is quantified to show how important it is to consider correcting the workfunction shifts (Upzc) along the reaction path.
5. Note the slope and beta in model 2c is potential-dependent. An effective beta can be calculated by averaging the beta over a potential range of interest


## Jupyter Notebook: Jupyter_Barriers_EDL.ipynb

--Current Capabilities--
1. Interactive GUI to input DFT data and model parameters for sensitivity analysis
2. Plots Free Energy Changes wrt. Potential (Potential - Potential of zero charge) given different dielectric constant and width of EDL in a Helmholtz Model
3. Analyzes the magnitude of different complexities of electrification on free energy changes given dielectric constant and width of EDL
4. Sensitivity Analysis due to approximate dielectric constants and widths of the EDL by model 2c via ipywidgets 
5. An interactive symmetry factor (beta) calculator and its sensitivity to approximated EDL widths EDL via ipywidgets
6. Optional gui window that exports data analysis as an excel sheet in the folder of script 

--Notes--
1. Currently the code isn't written to save dictionary results. I am working on to record previous entries via a dropdown menu. 
2. Beta Calculator is calculated for Cation+ transfer (H+) from bulk to the surface. Thus, beta builds from 0 to 1 as the EDL effects and changes in dipole moment/polarizability become more significant (1 is where dipole moment changes, polarizability changes, and EDL is not significant). Simply change e = 0 if you are studying a reaction where H+ is not needed with the transfer of an electron. Beta will then decrease from 1 to 0 as the EDL effects and changes in dipole moment/polarizability become more significant. 
3. Both the beta calculator and the sensitivity analysis ranges for the dielectric constant and the EDL width can be altered in the function. 
4. The excel sheet is exported to the same folder. I plan to add a feature that it has the option to make entries in different sheets and such.


## Bash Scripts

### polar

This is a simple script I made that sets up polarizability calculations for you given a directory with an optimized geometry. Simply, polar copys the VASP inputs (renames CONTCAR to POSCAR) and creates a directory called "field" where a set of singlepoint vasp calculations of efield from 0.1 to 0.6. These files are copied into each directory and the submission script is executed.
A few notes:
1. You need to change the submission script to match yours. Mine is called "SLURM.VASP". 
2. Feel free to change the script and iterate through the desired range of efield values. I find generally 0 to 0.5 with the five additional single points should suffice in determining the parabola of energy change in an electric field.
    a. I plan to create another script using gnuplot or python that can graph and process the value of the polarizability

### getenergies

This is a script generally nice to have when dealing with multiple VASP calculations in a directory. T
his script is executed as "getenergies path", where the path is the directory (and subdirectories) of interest. It 1) checks if the VASP jobs in the subdirectory are converged and 2) if so, grep the energy from the OUTCAR. 

## Additional tools

### QVASP and VASPKIT
For calculation the Workfunction (potential of zero charges), I recommend both QVASP and VASPKIT as they have an easy way to determine the WF ([Link here](https://sourceforge.net/projects/qvasp/)). If you want to analyze the xy average potenital w.r.t z of your surface, VASPKIT is built into QVASP and you can analyze this. 