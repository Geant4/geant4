import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data from Species.txt
fileName1 = "Species.txt"
data = pd.read_csv(fileName1, delim_whitespace=True, header=None, names=["Time", "Value", "Err", "Species"],
                   comment="#")

# Load data from e_aq.txt
fileName2 = "e_aq.txt"
data_eaq = pd.read_csv(fileName2, delim_whitespace=True, header=None, names=["Time", "Value", "Err"], comment="#")
data_eaq["Species"] = "e_aq^-1"  # Assign species name

# Load data from OH.txt
fileName3 = "OH.txt"
data_oh = pd.read_csv(fileName3, delim_whitespace=True, header=None, names=["Time", "Value", "Err"], comment="#")
data_oh["Species"] = "°OH^0"  # Assign species name

# Log scale transformation (avoid log(0) error)
data["Time"] = np.log10(data["Time"].replace(0, np.nan))
data_eaq["Time"] = np.log10(data_eaq["Time"].replace(0, np.nan))
data_oh["Time"] = np.log10(data_oh["Time"].replace(0, np.nan))

# Define species to plot with LaTeX formatting
species_list = ["°OH^0", "e_aq^-1", "H3O^1", "H2O2^0", "H_2^0", "H^0"]
labels = [
    r"$\mathrm{\cdot OH}$",  # OH radical
    r"$\mathrm{e^-_{aq}}$",  # Aqueous electron
    r"$\mathrm{H_3O^+}$",  # Hydronium ion
    r"$\mathrm{H_2O_2}$",  # Hydrogen peroxide
    r"$\mathrm{H_2}$",  # Molecular hydrogen
    r"$\mathrm{H^\bullet}$"  # Hydrogen radical
]
marker_styles = ["s", "^"]  # Square (e_aq), Triangle (OH)

# Create figure and axes
fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
axes = axes.flatten()

for i, (species, label) in enumerate(zip(species_list, labels)):
    ax = axes[i]
    subset = data[data["Species"] == species]

    if not subset.empty:
        ax.plot(subset["Time"], subset["Value"], linestyle='-', marker='o', label="Simulation", color='black')
        ax.errorbar(subset["Time"], subset["Value"], yerr=subset["Err"], fmt='o', capsize=5, color='black')

    # Add e_aq.txt as additional points
    if species == "e_aq^-1":
        ax.scatter(data_eaq["Time"], data_eaq["Value"], marker=marker_styles[0], color='red',
                   label=r"$\mathrm{e^-_{aq} \ exp}$")

    # Add OH.txt as additional points
    if species == "°OH^0":
        ax.scatter(data_oh["Time"], data_oh["Value"], marker=marker_styles[1], color='blue',
                   label=r"$\mathrm{\cdot OH \ exp}$")

    ax.set_title(label, fontsize=14)
    ax.set_yscale("linear")
    ax.grid(True)
    ax.legend()

# Set common labels
fig.text(0.5, 0.001, r"$\mathrm{Time \ (log(ps))}$", ha='center', fontsize=14)
fig.text(0.005, 0.5, r"$\mathrm{G(Species/100 \ eV)}$", va='center', rotation='vertical', fontsize=14)

plt.tight_layout()
plt.show()
