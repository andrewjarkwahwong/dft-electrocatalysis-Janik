--Background--
This part of the repository is to provide a tool to calculate potential dependent electrochemical activation barriers with electrochemical double layer (EDL) consideration. 
The main advantage of our approach is to test the sensitivity of DFT prediction activation barriers to different models that approximate the EDL.

I have provided several different tools within the repository and links based on your preference, where all tools compute similar sensitivity analysis.

Please read: (INSERT NAVEEN'S PAPER) for more details on the theory and derivation of our approach.

--Available Tools--
Jupyter Notebook: Barrier_Calculator_Reduction.ipyn for reduction and Barrier_Calculator_Oxidation.ipyn for oxidation
Google Colab: (INSERT LINK HERE) for the same notebooks listed above but in Google Colab 
Google Sheet: (INSERT REDUCTION LINK) and (INSERT OXIDATION LINK)


--Current Capabilities--
1. Interactive GUI to input DFT data and model parameters for sensitivity analysis
2. Plots Free Energy Changes wrt. Potential (Potential - Potential of zero charge) given different EDL models, dielectric constant, and width of EDL
3. Analyzes the magnitude of different complexities of electrification on free energy changes given dielectric constant and width of EDL
4. Sensitivity Analysis due to approximate dielectric constants and widths of the EDL via ipywidgets 
5. An interactive symmetry factor (beta) calculator and its sensitivity to approximated EDL widths EDL via ipywidgets
6. Optional gui window that exports data analysis as an excel sheet in the folder of script 


--Inputs--

Initial/Final States:
1. DFT Free Energies: Output of DFT energies with corrections (ex: ZPVE, TS) to Free Energies
2. Dipole Moment: Dipole Moments calculated from VASP (IDIPOL =3, LDIPOL = .TRUE.) for initial and final states
3. Polarizability: Refer to SI Figure # of (NAVEEN PAPER) of how to calculate the polarizability of the surface and the adsorbate

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
11. Relative Permittivity: Dielectric Constant of the media within the interface (er = 1 for vacuum or er = 78.4 for bulk water)
12. Width of the EDL: Width is defined as the distance between the electrode surface and the countercharge ions in angstrom 
13. Vacuum to NHE: The voltage correction from absolute scale to NHE scale (commonly as 4.2 V to 4.8 V but test for your system)
14. Solvation Free Energy Change (eV): This term incorporates the total solvation free energy change along the reaction path