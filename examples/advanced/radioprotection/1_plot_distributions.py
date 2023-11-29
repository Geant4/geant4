import numpy as np
import matplotlib.pyplot as plt
import glob

def Normalize(f, x) :
	dx = x[1:] - x[:-1]
	integral = np.sum(f[:-1] * dx)
	f_norm = f / integral
	return f_norm

# When running in MT and outputting to cvs, each thread outputs to a different file.
# Therefore I need to cycle through each thread and append them all together.

energy = np.empty([0,1])
length = np.empty([0,1])

ntuple_name = "radioprotection_nt_102"
output_name = ntuple_name + "_t" + "*" + ".csv"
output_list = glob.glob(output_name)

for this_file in output_list :	
	energy_thread, length_thread = np.loadtxt(this_file, delimiter=',', unpack=True, usecols=(0,1))
	
	energy = np.append(energy, energy_thread)
	length = np.append(length, length_thread)

# Experimentally, the mean path length is calculated geometrically as a mean chord length.
# Here for convenience it's taken by averaging the effective path lengths.
mean_path_length = np.average(length)

# A conversion factor can be used to convert the target material to water or tissue equivalent.
# Its choice depends on the material and can be calculated in different ways.
# It is suggested that the user replace the following value with his own.
conversion_factor = 1.

# If the previous value is not overriden by the user, this script will attempt to read geometry.mac
# and provide a factor accordingly
if conversion_factor == 1. :
	detector = np.array(["Diamond", "MicroDiamond", "Silicon", "SiliconBridge"])
	factor = np.array([0.32, 0.32, 0.57, 0.57])	# conversion factor based on material stopping power
	
	with open("geometry.mac") as search:
		for line in search:
			line = line.rstrip()  # remove new line
			
			for this in detector :
				match = "/geometrySetup/selectDetector " + this
				
				if match == line:
					conversion_factor = factor[ detector == this ][0]
				else:
					conversion_factor = factor[ detector == "Diamond" ][0]	# default detector type (no macro used)


y = energy * conversion_factor / mean_path_length


# The spectrum is now binned logarithmically, to avoid oscillations at higher energies
# (due to fewer counts) that wouldn't much meaning.

minimum = np.amin(y)
maximum = np.amax(y)

exp_start = np.floor(np.log10(minimum))
exp_end = np.ceil(np.log10(maximum))

n_decades = int(exp_end - exp_start)

# Number of logarithmic bins per decade:
# Higher values give better resolution, but lead to oscillations
# (especially at high energy) if your statistic has too few counts.
bins_per_dec = 60

n_bins = n_decades * bins_per_dec

y_bins = np.zeros(n_bins)

y_bins[0] = 10**exp_start

for i in range(1, n_bins) :
	y_bins[i] = y_bins[i-1] * 10**( 1 / bins_per_dec )


# Create the histogram

# For now f is a number of counts...
f = np.histogram( y, bins=y_bins ) [0]

tot_counts = np.sum(f)

# ... so now I turn f into a density
bin_width = y_bins[1:] - y_bins[:-1]
f = f / bin_width

f = np.append(f, 0.)	# give f and y_bins arrays the same size

# Normalize the spectra to unit area under the curve
f = Normalize(f, y_bins)
d = Normalize(y_bins*f, y_bins)


# Save to file

output_file = "analysed_spectra.csv"
header = "y[keV/um], f(y)[um/keV], d(y)[um/keV]"

np.savetxt( output_file, np.c_[ y_bins, f, d ], header=header, delimiter=',' )


# Plot

fig, ax1 = plt.subplots()

color = 'tab:blue'

ax1.semilogx( y_bins, y_bins * f, linewidth=0.5, color=color  )

ax1.set_xlabel(r'$y \,\, [keV / \mu m]$')
ax1.set_ylabel(r'$y \cdot f(y) $', color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'

ax2.semilogx( y_bins, y_bins * d, linewidth=0.5, color=color )

ax2.set_ylabel(r'$y \cdot d(y) $', color=color)
ax2.tick_params(axis='y', labelcolor=color)

title = str(tot_counts) + " counts, " + str(bins_per_dec) + " bins per decade, " + str(conversion_factor) + " conversion factor"
fig.suptitle(title)

fig.tight_layout()

plt.subplots_adjust(top=0.92)
plt.show()
