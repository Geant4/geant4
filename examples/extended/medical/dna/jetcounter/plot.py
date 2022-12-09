#!/usr/bin/env python3
# coding: utf-8
# *********************************************************************
# To execute this script simply type at your machine's prompt:
#   python plot.py
# OR
#   python3 plot.py
# This script needs:
# 1. A bunch of CShistory_t* text files created during simulation.
# 2. (optionally) energy.spectrum file provided with the example.
# *********************************************************************

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')  # to work also in sessions without display
import matplotlib.pyplot as plt


fs = []
for filename in os.listdir():
    if "CShistory_t" in filename:
        fs.append(np.loadtxt(filename))
events_data = np.concatenate(fs)

nevents = events_data.shape[0]

# ###
# ICSD plot

cluster_sizes = events_data[:, 0]

M1 = cluster_sizes.mean()

bin_edges = np.arange(cluster_sizes.max() + 2)

y, bin_edges = np.histogram(cluster_sizes, bins=bin_edges, density=True)

icsd_fig = plt.figure(figsize=(10, 5))
plt.errorbar(bin_edges[:-1], y, yerr=np.sqrt(y / nevents),
             fmt='d', label=f"M₁ = {M1:.3g}")

plt.ylim(.3/cluster_sizes.size, 1)
plt.yscale("log")

plt.title("ICSD")
plt.xlabel("cluster size")
plt.ylabel("probability")
plt.legend()

icsd_fig.savefig("ICSD_py.png")


# ###
# energy spectra
initial_energies = events_data[:, 1]
interaction_energies = events_data[:, 2]
final_energies = events_data[:, 3]

labels = ["initial", "interaction", "final"]

spec_fig = plt.figure(figsize=(10, 5))
for i, energies in enumerate((initial_energies, interaction_energies, final_energies)):
    y, bin_edges = np.histogram(energies, bins=int(np.sqrt(nevents)), density=True)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    plt.errorbar(bin_centers, y, yerr=np.sqrt(y / nevents), drawstyle='steps-mid',
                 label="{} <E> = {:.3g} MeV".format(labels[i], energies.mean()))


try:
    source_spectrum_data = np.loadtxt("energy.spectrum")
    length, gain, offset = source_spectrum_data[:3]
    cumulative_counts = source_spectrum_data[3:]
    counts = np.diff(cumulative_counts)
    bins = np.arange(length-1)*gain+offset

    normalization_factor = counts.sum()*gain
    density = counts/normalization_factor

    mean_energy = np.sum(bins*density)*gain

    plt.plot(bins, density, drawstyle='steps-mid', label="input <E> = {:.3g} MeV".format(mean_energy))
except FileNotFoundError:
    print("Input spectrum file 'energy.spectrum' not found!")
    pass


plt.title("Energy spectra")
plt.xlabel("energy [MeV]")
plt.ylabel("probability density [1/MeV]")
plt.legend()

spec_fig.savefig("energy_spectra.png")
