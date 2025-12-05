#!/usr/bin/env python
# coding: utf-8

# Read and plot the simulation results of particle interactions in Oriented Crystals
# obatined through example ch2, which is baed on G4ChannelingFastSimModel.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import uproot
from matplotlib.colors import LogNorm  # optional, for log color scaling

#################################### INPUT FILE #########################################
# Set path and filename of the simulation file
G4_sim_path = ""
root_file = "results"

# Set whether to save plots (without displaying them) or just display them
save_fig = True
fig_path = G4_sim_path
#########################################################################################

# Create directory where to strore the figures if it does not exist
if fig_path != '' and not os.path.exists(fig_path):
    os.makedirs(fig_path)
    print('created fig_path:', fig_path)

# Open the simulation output root file
rf = uproot.open(G4_sim_path + root_file + '.root')
rf_content = [item.split(';')[0] for item in rf.keys()]
print('rf_content:', rf_content, '\n')

# Import the scoring ntuples and convert them into pandas dataframes
branches = ["eventID", "volume", "x", "y", "angle_x", "angle_y", \
            "Ekin" , "particle", "particleID", "parentID"]
branchesprimary = branches + ["incoming_angle_x", "deflection_angle_x", \
                              "incoming_angle_y", "deflection_angle_y"]

df_in = rf['crystal'].arrays(branches, library='pd')
df_prim = rf['detector_primaries'].arrays(branchesprimary, library='pd')
df_ph = rf['detector_photons'].arrays(branches, library='pd')
df_sec = rf['detector_secondaries'].arrays(branches, library='pd')
df_missed = rf['missed_crystal'].arrays(branches, library='pd')

#########################################################################################
# Plot angle_x distribution of primaries at the detector after interaction with a crystal

############# INPUT #############
# Feel free to modify according to your needs
Nmax = 100000000 #max number of events to elaborate

# Feel free to replace df_prim by df_in, df_ph, df_sec or df_missed
# Feel free to replace "angle_x" by other ntuples from branches and
# from branchesprimary (for df_prim)
# ONLY NUMERIC VALUES
datax = df_prim["angle_x"][:Nmax]*1.e3  #mrad <= rad (feel free to modify the coefficient)

# Feel free to modify the number of bins and the plot range
NbinTheta = 100
rangeTheta = [-1, 2] #mrad

# Set whether to use linear o log scale
use_log_y = False  # set True for LogNorm color scale

# Feel free to modify the names of axes
plt_xlabel = '$\\theta_x$ [mrad]'
plt_ylabel = 'PDF:  1/N dN/d$\\theta_x$ [mrad]$^{-1}$'

# Feel free to modify the filename to save the plot
filename = 'thetaXdistribution.pdf'

#some plt parameters
fs = 16
lw = 2
#################################

# Create 1D histogram
thetaXdistrib, thetaEdges = np.histogram(datax.values, \
                                         bins=NbinTheta, range=rangeTheta, density=True)
thetabin = thetaEdges[:-1] + (thetaEdges[1]-thetaEdges[0])*0.5
plt.figure(figsize=(9, 6))
plt.grid()
plt.plot(thetabin, thetaXdistrib, linewidth=lw, alpha=1, label='')
plt.xlabel(plt_xlabel, fontsize=fs)
plt.ylabel(plt_ylabel, fontsize=fs)

# Set log scale
if use_log_y:
    plt.yscale('log',base=2)

# Save the plot or just show it
if save_fig:
    plt.savefig(fig_path + filename)
    plt.close()

#########################################################################################
# angle_x_in - angle_x_defl distribution of primaries at the detector after interaction with a crystal

############# INPUT #############
# Feel free to modify according to your needs
Nmax = 100000000 #max number of events to elaborate

# Example data (replace these with your real arrays)
# datatetaxin and datatetadeflx must be the same length
datatetaxin = df_prim["incoming_angle_x"][:Nmax]*1.e3  #mrad <= rad (feel free to modify the coefficient)
datatetadeflx = df_prim["deflection_angle_x"][:Nmax]*1.e3  #mrad <= rad (feel free to modify the coefficient)

# Feel free to modify the number of bins and the plot range
NbinTheta = 50

# Feel free to modify the plot range
xrange = (-0.1, 0.1)
yrange = (-1, 2)

# Set whether to use linear o log scale
use_log_color = True  # set True for LogNorm color scale

# Feel free to modify the names of axes
plt_xlabel2 = '$\\theta_{x in}$ [mrad]'
plt_ylabel2 = '$\\theta_{x defl}$ [mrad]'

# Feel free to modify the filename to save the plot
filename2 = 'thetaXin_thetaXdefl.pdf'

#some plt parameters
fs = 16
lw = 2
#################################

# Create 2D histogram
plt.figure(figsize=(8, 6))
hist = plt.hist2d(
    datatetaxin,
    datatetadeflx,
    bins=NbinTheta,
    density=True,
    range=[xrange, yrange],
    norm=LogNorm() if use_log_color else None,
    cmap='jet'
)

# Add colorbar (PDF scale)
cbar = plt.colorbar()
cbar.set_label('PDF', fontsize=fs)

# Labels and title
plt.xlabel(plt_xlabel2, fontsize=fs)
plt.ylabel(plt_ylabel2, fontsize=fs)

# Save the plot or just show it
if save_fig:
    plt.savefig(fig_path + filename2)
    plt.close()
        
################################################################################################
# Plot spectrum

############# INPUT #############
# Feel free to modify the collimator parameters
apply_collimation = True
coll_angle = 2.3183 #mrad

# Feel free to modify
NbinE = 20
rangeE = [0, 10] #MeV

# path of the spectrum file obtained using all the Bair-Katkov integration photons
BK_spectrum_file = "Spectrum.dat"
#################################

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
spectrum0, EbinEdges = np.histogram(Eph, bins=NbinE, range=rangeE, density=False)
Ebin = EbinEdges[:-1] + (EbinEdges[1]-EbinEdges[0])*0.5
stepx = Ebin[1]-Ebin[0]
Nprimaries = df_in["Ekin"].size

spectrum = spectrum0 / (Nprimaries*stepx)
spectral_intensity = Ebin * spectrum

# Statistical uncertainties: sqrt(N)
spectrum_err = np.sqrt(spectrum0) / (Nprimaries * stepx)
spectral_intensity_err = Ebin * spectrum_err

# Read the spectrum file obtained using all the Bair-Katkov integration photons
BK_spectrum = np.loadtxt(G4_sim_path+BK_spectrum_file, dtype='float', comments='#', \
                         delimiter=' ', skiprows=1, unpack=True)
E_ext = BK_spectrum[0]
S_ext = BK_spectrum[1]

# Plot the photon energy spectrum
fig = plt.figure(figsize=(13, 6))
fs = 16
lw = 2
bw = 0.6
color0 = '#B1B3FB'

plt.subplot(1,2,1)
plt.bar(Ebin, spectrum, width=bw, color=color0, linewidth=lw, alpha=1, label='secondary photons')
plt.errorbar(Ebin, spectrum, yerr=spectrum_err, fmt='o', color='k', capsize=3)
plt.xlim(rangeE)
plt.plot(E_ext, S_ext, 'r-', lw=2.5, label='from '+BK_spectrum_file)
plt.title('Emitted photon spectrum')
plt.xlabel('E [MeV]', fontsize=fs)
plt.ylabel('1/N dN/dE', fontsize=fs)
plt.legend()
#plt.yscale('log')

plt.subplot(1,2,2)
plt.bar(Ebin, spectral_intensity, width=bw, color=color0, linewidth=lw, alpha=1, label='secondary photons')
plt.errorbar(Ebin, spectral_intensity, yerr=spectral_intensity_err, fmt='o', color='k', capsize=3)
plt.plot(E_ext, E_ext * S_ext, 'r-', lw=2.5, label='from '+BK_spectrum_file)
plt.title('Emitted photon spectral intensity')
plt.xlabel('E [MeV]', fontsize=fs)
plt.ylabel('1/N dW/dE', fontsize=fs)
plt.xlim(rangeE)
plt.legend()
#plt.yscale('log')
if save_fig:
    plt.savefig(fig_path + 'spectrum.pdf')
    plt.close()
