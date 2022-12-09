import numpy as np
from numpy import exp, sqrt, log, pi

def Integrate(f, x) :
	dx = x[1:] - x[:-1]
	integral = np.sum(f[:-1] * dx)
	return integral

input_file = "analysed_spectra.csv"
y, f, d = np.loadtxt(input_file, delimiter=',', unpack=True)


# --------------------------------------
#        Microdosimetric means
# --------------------------------------

yBarF = Integrate(y*f, y)
yBarD = Integrate(y*d, y)

print("")
print("Weighted averages of:")
print("f(y):   ", round(yBarF, 3) )
print("d(y):   ", round(yBarD, 3) )
print("")


# --------------------------------------
#        Quality factor
# --------------------------------------

# Quality factor of the radiation, used in radiation protection to weight the absorbed dose.

# The approximation employed here follows ICRU 40.
Q_y = 5510./y * ( 1 - exp(-5E-5*y**2 -2E-7*y**3) )

Q = Integrate(Q_y*d, y)

print("Average quality factor Q:   ", round(Q, 3) )
print("")


# --------------------------------------
# RBE via Microdosimetric Kinetic Model
# --------------------------------------

# The MKM can calculate the alpha parameter of the survival curve from an average of the microdosimetric
# spectrum, weighted by y and corrected for saturation.

# Saturation coefficient:
y0 = 125. # keV/um	according to ICRU 36
#y0 = 150. # keV/um	according to Kase 2006, for Carbon ions

# Saturation corrected average, usually written as y*:
y_star = y0**2 * Integrate( (1 - exp(-(y/y0)**2))* f, y ) / yBarF

print("y*:   ", round(y_star, 3) )

# The following parameters depend on the cell line, and need to be known in advance.
# The ones below are the ones employed by Kase 2006 for HSG cells.
# Refer to said paper for a description of each parameter.
alpha_0 = 0.13 # Gy**-1
beta = 0.05 # Gy**-2
rho = 1. # g/cm**3
r_d = 0.42 # um

keV_to_Gy=1.602E-16; um_to_cm=1.E-4; g_to_kg=1.E-3
rho*=g_to_kg; r_d*=um_to_cm; y_star*=(keV_to_Gy/um_to_cm)

alpha = alpha_0 + beta/(rho*pi*r_d**2) * y_star

print("alpha:   ", round(alpha, 3) )

# The RBE is not an absolute value, but depends on the desired survival fraction S
# and on the corresponding photon dose. The latter is also taken from Kase 2006.
S = 0.1
Dx = 5 # Gy	for S=0.1

RBE_MKM = Dx*2*beta / ( sqrt(alpha**2-4*beta*log(S)) - alpha )

print("RBE via MKM (HSG cells):   ", round(RBE_MKM, 3) )
print("")


# --------------------------------------
#   RBE via Loncol weight function r(y)
# --------------------------------------

# A weight function can be used to estimate the RBE for S=0.1 by weighting the whole spectrum.
# The one by Loncol 1994 provides a good estimate for neutron or proton beams targeting crypt cells,
# but not necessarily for other beams-target combinations.

weight_function = "weight_function.csv"
y_wf, r_wf, stdev_r_wf = np.loadtxt(weight_function, delimiter=',', unpack=True)

# The standard deviation of r(y) is due to the uncertainty of radiobiological measurements,
# and as such doesn't represent a limitation of the model.
# Consequently it won't be used in the following RBE estimate.

# The abscissae of f(y) and r(y) don't (usually) match
# For every value of f(y), I find the corresponding r(y) via linear interpolation
size = len(y)
r = np.zeros(size)

for i in range(size) :
	j = np.argmax( y_wf > y[i] )	# index of first element of y_wf greater than y[i]
	
	r[i] = ( r_wf[j] + r_wf[j-1] ) /2

RBE_Lwf = Integrate(r*d, y)

print("RBE via weight function (crypt cells):   ", round(RBE_Lwf, 3) )
print("")
