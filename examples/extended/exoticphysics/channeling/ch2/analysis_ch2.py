#!/usr/bin/env python
# coding: utf-8

# Read and plot the simulation results of particle interactions in Oriented Crystals
# obatined through example ch2, which is baed on G4ChannelingFastSimModel.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import uproot

################################### INPUT ############################################
# Set path and filename of the simulation file
G4_sim_path = ""
root_file = "results"

Nmax = 1e5 #max number of events to elaborate

save_fig = True
fig_path = G4_sim_path

apply_collimation = False
coll_angle = 20/8627 #rad
NbinE = 25
rangeE = [0, 10] #MeV

NbinTheta = 100
rangeTheta = [-1, 2] #mrad
######################################################################################

# Create figure directory if it does not exist
if fig_path != '' and not os.path.exists(fig_path):
    os.makedirs(fig_path)
    print('created fig_path:', fig_path)

# Open the simulation output root file
rf = uproot.open(G4_sim_path + root_file + '.root')
rf_content = [item.split(';')[0] for item in rf.keys()]
print('rf_content:', rf_content, '\n')

# Import the scoring ntuples and convert them into pandas dataframes
branches = ["eventID", "volume", "x", "y", "angle_x", "angle_y", 
            "Ekin" , "particle", "particleID", "parentID"]
df_in = rf['crystal'].arrays(branches, library='pd')
df_out = rf['detector'].arrays(branches, library='pd')
df_ph = rf['detector_photons'].arrays(branches, library='pd')

# Define in and out dataframes 
df_in_all_primary = df_in[df_in.parentID == 0]
df_out_all_primary = df_out[df_out.parentID == 0]
Nmax = min([int(Nmax), len(df_out_all_primary)])
df_in_primary = df_in_all_primary[:Nmax]
df_out_primary = df_out_all_primary[:Nmax]

# Select only the columns useful for deflection
df_in_primary_sel = df_in_primary[["eventID", "angle_x", "angle_y"]]
df_out_primary_sel = df_out_primary[["eventID", "angle_x", "angle_y"]]
del df_in_primary, df_out_primary

# Array with photon energies and angles
Eph = df_ph['Ekin'].values #MeV 
Nph = len(Eph)
print("number of emitted photons:", Nph)
thetaX_ph = df_ph['angle_x'].values*1e3 #rad -> mrad
thetaY_ph = df_ph['angle_y'].values*1e3 #rad -> mrad

# Take only the photons inside the collimator acceptance
theta_ph = np.sqrt(thetaX_ph**2 + thetaY_ph**2)   
if apply_collimation: 
    thetaX_ph = thetaX_ph[theta_ph <= coll_angle]
    thetaY_ph = thetaY_ph[theta_ph <= coll_angle]
    Eph = Eph[theta_ph <= coll_angle]
    theta_ph = theta_ph[theta_ph <= coll_angle]
    
# Calculate the scored photon energy spectrum
spectrum, EbinEdges = np.histogram(Eph, bins=NbinE, range=rangeE, density=True)
Ebin = EbinEdges[:-1] + (EbinEdges[1]-EbinEdges[0])*0.5
spectral_intensity = Ebin * spectrum

# Plot the photon energy spectrum
fig = plt.figure(figsize=(13, 6))
fs = 16
lw = 2
bw = 0.6
plt.subplot(1,2,1)
plt.bar(Ebin, spectrum, width=bw, linewidth=lw, alpha=1, label='')
plt.title('Emitted photon spectrum')
plt.xlabel('E [MeV]', fontsize=fs)
plt.ylabel('1/N$\\times$dN/dE', fontsize=fs)
plt.yscale('log')
plt.subplot(1,2,2)
plt.bar(Ebin, spectral_intensity, width=bw, linewidth=lw, alpha=1, label='')
plt.title('Emitted photon spectral intensity')
plt.xlabel('E [MeV]', fontsize=fs)
plt.ylabel('1/N$\\times$dW/dE', fontsize=fs)
plt.yscale('log')
if save_fig:
    plt.savefig(fig_path + 'spectrum.jpg')
plt.close() 

# Plot angle_x distribution at the detector
thetaXdistrib, thetaEdges = np.histogram(df_out_primary_sel["angle_x"].values*1e3, \
                                         bins=NbinTheta, range=rangeTheta, density=True)
thetabin = thetaEdges[:-1] + (thetaEdges[1]-thetaEdges[0])*0.5
plt.figure(figsize=(9, 6))
plt.plot(thetabin, thetaXdistrib, linewidth=lw, alpha=1, label='')
plt.xlabel('$\\theta_X$ [mrad]', fontsize=fs)
plt.ylabel('1/N$\\times$dN/d$\\theta_X$', fontsize=fs)
if save_fig:
    plt.savefig(fig_path + 'thetaXdistribution.jpg')
plt.close() 

