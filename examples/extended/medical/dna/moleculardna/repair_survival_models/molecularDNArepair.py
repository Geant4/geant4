#!/usr/bin/env python
# coding: utf-8

# To run, change the links and definitions in this script and 
# type in the terminal containing the file 
# $python3 molecularDNArepair.py

# Authors: 
# K. Chatzipapas (*), D. Sakata, H. N. Tran, M. Dordevic, S. Incerti et al.
# (*) contact: k.chatzipapas@yahoo.com

# Define filenames
outputFile = "molecularDNArepair.txt"
iRootFile  = "../molecular-dna.root"
print("\nInput file:",iRootFile,"\n")

# Define cell parameters
r3 = 7100*1e-09 * 2500*1e-09 * 7100*1e-09 # a * b * c   / Chromosome size, as defined in the mac file, but in meters. If sphere, a=b=c
NBP = 6405886128 # bp / Length of the DNA chain in Mbp
mass = 997 * 4 * 3.141592 * r3 / 3     # waterDensity * 4/3 * pi * r3 in kg

# Integration parameters
t0 = 0.0;   # Starting time point (dimensionless)
t1 = 25;  # Final time point (dimensionless)
dt = 0.001 #0.001; # Intergration time step (dimensionless)
numpoints = int( (t1 - t0) / dt )


# ONLY ADVANCED USERS AFTER THIS POINT
# Repair model parameters may be find at line 105-172
#####################################################################################################################
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint as od

# Initialization
DSBBPID = []
i = 0
repDSB = 0
irrDSB = 0
acc_edep = 0
dose = 0
eVtoJ = 1.60218e-19


f = root.TFile(iRootFile)

fTree = f.Get("tuples/damage")
for e in fTree:
    if (e.TypeClassification=="DSB" or e.TypeClassification=="DSB+" or e.TypeClassification=="DSB++"):
        DSBBPID.append((e.Event, e.BasePair))
        i += 1

gTree = f.Get("tuples/classification")
for e in gTree:
    repDSB += e.DSB
    irrDSB += e.DSBp + e.DSBpp  # 2*eDSBpp   could also be used
    
ffTree = f.Get("tuples/chromosome_hits")
for ev in ffTree:
    acc_edep += (ev.e_chromosome_kev + ev.e_dna_kev) *1e3  # eV    



# Calculate the total number of DSBs
totalDSB = repDSB + irrDSB

# Calculate the absorbed dose
dose = acc_edep * eVtoJ / mass    # Dose in Gy

#print ("r3=", r3)
#print("Dose :", dose, "\ni :", i, "\nTotal DSB :", totalDSB, "\nRepairable DSBs :", repDSB, "\nIrrepairable DSBs :", irrDSB)

NBP_of_cell = 5000000000   # This is a constant of the model, as defined by Belov et al.
NBPcell = NBP/NBP_of_cell
totalDSBperGyCell = totalDSB/NBPcell/dose
print ("Dose :", dose, "\nTotal DSB/Gy/Cell :", totalDSBperGyCell)
S1 = repDSB/NBPcell/dose
S2 = irrDSB/NBPcell/dose
print("Repairable DSBs/Gy/cell :", S1, "\nIrrepairable DSBs/Gy/cell :", S2)


# Model parameters
fDz     = 1 # Dose is set to 1 Gy, because falpha is normalized to DSB/Gy/Cell
falpha  = totalDSBperGyCell #27.5 #totalDSB # total DSB
fNirrep = irrDSB/totalDSB # irrepDSB # racirreparable = (G4double)fIrreparableDSBs/(G4double)fTotDSBs;
if (fNirrep == 0):
    fNirrep = 0.01   # For not printing error
print("Ratio of irrepairable to repairable damage :", fNirrep)

# Initial conditions for ODE
NbEquat = 29
Y  = [0] * NbEquat
#print (Y[1])
Y[0] = falpha * fDz
YP = [0] * NbEquat
print("\nY[0] used in the model :", Y[0])

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t = [(t1-t0) * float(i) / (numpoints - 1) for i in range(numpoints)]


#  DIMENSIONAL REACTION RATES DEFINITON (check Belov et al.) (Belov OV, et al. J Theor Biol. 2015 Feb 7;366:115-30. doi: 10.1016/j.jtbi.2014.09.024)
#
#------------NHEJ--------------
K1 = 11.052 #other:10.02         #11.052             # M-1*h-1  # change parametrs based on the article of Belov et al.
Kmin1 = 6.6*1e-01       #6.59999*1e-04   # h-1   # other 6.6*1e-01       #
K2 = 5.82*1e+05   #18.8305*(1.08517-np.exp(-21.418/pow(fDz,1.822))) # M-1*h-1  # other:5.82*1e+05    #
#print (K2)
Kmin2 = 5.26*1e-01      #h-1
K3 = 1.86               # h-1
K4 = 1.38*1e+06         # M-1*h-1
Kmin4 = 3.86*1e-04      # h-1
K5 = 15.24              # M-1*h-1
Kmin5 = 8.28            # h-1
K6 = 18.06              # M-1*h-1
Kmin6 = 1.33            # h-1
K7 = 2.73*1e+05         # M-1*h-1
Kmin7 = 3.2             # h-1
#K8 is a normalization factor and sometimes it needs to be adjusted.
normK8 = 0.35
K8 = normK8*5.52*1e-01  # h-1
K9 = 1.66*1e-01         # h-1
K10 = (1.93*1e-07)/fNirrep # M
K11 = 7.50*1e-02        # h-1
K12 = 11.1              # h-1
#
#------------HR--------------
P1 = 1.75*1e+03         # M-1*h-1
Pmin1 = 1.33*1e-04      # h-1
P2 = 0.39192            # h-1
Pmin2 = 2.7605512*1e+02 # h-1
P3 = 1.37*1e+04         # M-1*h-1
Pmin3 = 2.34            # h-1
P4 = 5.52*1e-2     #3.588*1e-02        # h-1
P5 = 1.20*1e+05         # M-1*h-1
Pmin5 = 8.82*1e-05      # h-1
P6 = 1.87*1e+05    #1.54368*1e+06      # M-1*h-1
Pmin6 = 1.55*1e-03      # h-1
P7 = 21.36         #1.4904             # h-1
P8 = 1.20*1e+04         # M-1*h-1
Pmin8 = 2.49*1e-04      # h-1
P9 = 4.88*1e-1     #1.104              # h-1
P10 = 7.20*1e-03        # h-1
P11 = 6.06*1e-04        # h-1
P12 = 2.76*1e-01        # h-1
#
#------------SSA--------------
Q1 = 7.8*1e+03      #1.9941*1e+05       # M-1*h-1
Qmin1 = 1.71*1e-04      # h-1
Q2 = 3*1e+04        #4.8052*1e+04       # M-1*h-1
Q3 = 6*1e+03            # M-1*h-1
Qmin3 = 6.06*1e-04      # h-1
Q4 = 1.66*1e-06     #1.62*1e-03         # h-1
Q5 = 8.40*1e+04         # M-1*h-1
Qmin5 = 4.75*1e-04      # h-1
Q6 = 11.58              # h-1
#
#-------alt-NHEJ (MMEJ)--------
R1 = 2.39*1e+03         # M-1*h-1
Rmin1 = 12.63           # h-1
R2 = 4.07*1e+04         # M-1*h-1
R3 = 9.82               # h-1
R4 = 1.47*1e+05         # M-1*h-1
Rmin4 = 2.72            # h-1
R5 = 1.65*1e-01         # h-1
#
# Scalling rate XX1
XX1 = 9.19*1e-07 # M
XX3 = 2.3 *1e-12 # M
#
# DIMENSIONLESS REACTION RATES 
#
#------------NHEJ--------------
k1 = K1*XX1/K8
kmin1 = Kmin1/K8
k2 = K2*XX1/K8 
kmin2 = Kmin2/K8     
k3 = K3/K8               
k4 = K4*XX1/K8         
kmin4 = Kmin4/K8      
k5 = K5*XX1/K8             
kmin5 = Kmin5/K8            
k6 = K6*XX1/K8              
kmin6 = Kmin6/K8            
k7 = K7*XX1/K8          
kmin7 = Kmin7/K8           
k8 = K8/K8         
k9 = K9/K8         
k10 = K10/XX1 
k11 = K11/K8       
k12 = K12/K8            
#
#------------HR--------------
p1 = P1*XX1/K8       
pmin1 = Pmin1/K8     
p2 = P2/K8           
pmin2 = Pmin2/K8  
p3 = P3*XX1/K8          
pmin3 = Pmin3/K8        
p4 = P4/K8       
p5 = P5*XX1/K8       
pmin5 = Pmin5/K8       
p6 = P6*XX1/K8        
pmin6 = Pmin6/K8       
p7 = P7/K8              
p8 = P8*XX1/K8          
pmin8 = Pmin8/K8      
p9 = P9/K8          
p10 = P10/K8        
p11 = P11/K8        
p12 = P12/K8       
#
#------------SSA--------------
q1= Q1*XX1/K8       
qmin1 = Qmin1/K8      
q2 = Q2*XX1/K8       
q3 = Q3*XX1/K8        
qmin3 = Qmin3/K8     
q4 = Q4/K8         
q5 = Q5*XX1/K8          
qmin5 = Qmin5/K8     
q6 = Q6/K8             
#
#-------alt-NHEJ (MMEJ)--------
r1 = R1*XX1/K8         
rmin1 = Rmin1/K8        
r2 = R2*XX1/K8         
r3 = R3/K8            
r4 = R4*XX1/K8         
rmin4 = Rmin4/K8           
r5 = R5/K8        
#------------------------------------

p = [ k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12,
      kmin1, kmin2, kmin4, kmin5, kmin6, kmin7,
      p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,
      pmin1, pmin2, pmin3, pmin5, pmin6, pmin8,
      q1, q2, q3, q4, q5, q6,
      qmin1, qmin3, qmin5,
      r1, r2, r3, r4, r5,
      rmin1, rmin4
      ]


# Model definition
def DSBRepairPathways( Y , t , p ):
    # SYSTEM OF DIFFERENTIAL EQUATIONS
    
    # Concentrations of repair enzymes set to be constant
    # X1  # [Ku]      # X2  # [DNAPKcsArt]   # X3  # [LigIV/XRCC4/XLF]
    # X4  # [PNKP]    # X5  # [Pol]          # X6  # [H2AX]
    # X7  # [MRN/CtIP/ExoI/Dna2]             # X8  # [ATM]
    # X9  # [RPA]     # X10 # [Rad51/Rad51par/BRCA2]
    # X11 # [DNAinc]  # X12 # [Rad52]        # X13 # [ERCC1/XPF] 
    # X14 # [LigIII]  # X15 # [PARP1]        # X16 # [Pol]
    # X17 # [LigI]

    X1= X2= X3= X4= X5= X6= X7= X8= X9= X10= X11= X12= X13= X14= X15= X16= X17=  400000. 
    
    # ----- NHEJ ----------
    YP[0] = fNirrep - k1*Y[0]*X1 + kmin1*Y[1] - p1*Y[0]*X1 + pmin1*Y[10] # [DSB]
    YP[1] = k1*Y[0]*X1 - kmin1*Y[1] - k2*Y[1]*X2 + kmin2*Y[2] # [DBS * Ku]   
    YP[2] = k2*Y[1]*X2 - k3*Y[2] - kmin2*Y[2] # [DSB * DNA-PK/Art]
    YP[3] = k3*Y[2] - k4*(Y[3]*Y[3]) + kmin4*Y[4] # [DSB * DNA-PK/ArtP]
    YP[4] = k4*(Y[3]*Y[3]) - kmin4*Y[4] - k5*Y[4]*X3 + kmin5*Y[5] # [Bridge]
    YP[5] = kmin6*Y[6] + k5*Y[4]*X3 - kmin5*Y[5] - k6*Y[5]*X4 # [Bridge * LigIV/XRCC4/XLF]
    YP[6] = -kmin6*Y[6] - k7*Y[6]*X5 +  kmin7*Y[7] + k6*Y[5]*X4 # [Bridge * LigIV/XRCC4/XLF * PNKP]
    YP[7] = k7*Y[6]*X5 - k8*Y[7] - kmin7*Y[7] # [Bridge * LigIV/XRCC4/XLF * PNKP * Pol]
    YP[8] = r5*Y[28] + k8*Y[7] + p12*Y[18] + p11*Y[19] + q6*Y[24] # [dsDNA]
    YP[9] = (k9* (Y[3]+Y[4]+Y[5]+Y[6]+Y[7]) * X6)/ (k10+ Y[3]+ Y[4]+ Y[5]+ Y[6]+ Y[7]) - k11*Y[8]- k12*Y[9] # [gH2AX foci]

    # ----- HR ---------- 
    YP[10] = p1*Y[0]*X7 - pmin1*Y[10] - p3*Y[10]*Y[11] + pmin3*Y[12]     # [MRN/CtIP/ExoI/Dna2]
    YP[11] = p2*X8 - pmin2*Y[11] - p3*Y[10]*Y[11] + p4*Y[12] + pmin3*Y[12]     # [ATMP]
    YP[12] = p3*Y[10]*Y[11] - p4*Y[12] - pmin3*Y[12]     # [DSB * MRN/CtIP/ExoI/Dna2 * ATMP]
    YP[13] = rmin1*Y[25] + p4*Y[12] - r1*X15*Y[13] - p5*Y[13]*X9 + pmin5*Y[14]     # [ssDNA]
    YP[14] = pmin6*Y[15] + p5*Y[13]*X9 - pmin5*Y[14] - p6*Y[14]*X10 - q1*Y[14]*X12 + qmin1*Y[20]     # [ssDNA * RPA]
    YP[15] = -p7*Y[15] - pmin6*Y[15] + p6*Y[14]*X10    # [ssDNA * RPA * Rad51/Rad51par/BRCA2]
    YP[16] = p7*Y[15] - p8*Y[16]*X11 + pmin8*Y[17] # [Rad51 filament]
    YP[17] = p8*Y[16]*X11 - p9*Y[17] - pmin8*Y[17] # [Rad51 filament * DNAinc]
    YP[18] = p9*Y[17] - p10*Y[18] - p12*Y[18] # [D-loop]
    YP[19] = p10*Y[18] - p11*Y[19] # [dHJ]

    # ----- SSA ---------- 
    YP[20] = q1*Y[14]*X12 - qmin1*Y[20] -  q2*(Y[20]*Y[20])    # [ssDNA * RPA * Rad52]
    YP[21] = q2*(Y[20]*Y[20]) - q3*Y[21]*X13 + qmin3*Y[22] # [Flap]
    YP[22] = q3*Y[21]*X13 - q4*Y[22] - qmin3*Y[22] # [Flap * ERCC1/XPF]
    YP[23] = q4*Y[22] - q5*Y[23]*X14 + qmin5*Y[24] # [dsDNAnicks]
    YP[24] = q5*Y[23]*X14 - q6*Y[24] - qmin5*Y[24] # [dsDNAnicks * LigIII]

    # ----- MMEJ ---------- 
    YP[25] = -rmin1*Y[25] - r2*Y[25]*X16 + r1*X15*Y[13] # [ssDNA * PARP1]
    YP[26] = r2*Y[25]*X16 - r3*Y[26] # [ssDNA * Pol]
    YP[27] = r3*Y[26] - r4*Y[27]*X17 + rmin4*Y[28] # [MicroHomol]
    YP[28] = r4*Y[27]*X17 - r5*Y[28] - rmin4*Y[28] # [MicroHomol * LigI]

    #---------------------------------------------------

    return YP

# Model execution
steps = od( DSBRepairPathways, Y , t , args=(p,) )

print("\nCalculation concluded. Check the final images and close them to end the program.!")

# Plotting results
KU = []
DNAPKcs = []
RPA = []
Rad51 = []
gH2AX = []
gH2AXnorm = []
i = 0
for i in range(len(t)):
    KU.append(steps[i][1]) 
    DNAPKcs.append(steps[i][3]+steps[i][4]+steps[i][5]+steps[i][6]+steps[i][7]) 
    RPA.append(steps[i][14]+steps[i][15]+steps[i][20]) 
    Rad51.append(steps[i][15]+steps[i][16]+steps[i][17]) 
    gH2AX.append(steps[i][9]) 
    gH2AXnorm.append(100*steps[i][9]/repDSB) 


plt.figure("KU")
plt.plot(t, KU)
plt.show()


plt.figure("DNAPKcs")
plt.plot(t, DNAPKcs)
plt.show()


plt.figure("RPA")
plt.plot(t, RPA)
plt.show()


plt.figure("Rad51")
plt.plot(t, Rad51)
plt.show()


#plt.xlim(0,25)
plt.figure("gH2AX")
plt.plot(t, gH2AX)
plt.savefig("g-H2AX.pdf")

with open(outputFile, 'w') as ofile:
    for m in range(len(gH2AX)):
        ofile.write(str(t[m])+ "\t" + str(gH2AX[m])+"\n")  

plt.show()


# Normalization, and plotting with available data, as published in Chatzipapas et al., arxiv https://doi.org/10.48550/arXiv.2210.01564
#plt.xlim(0,25)
plt.figure("gH2AXnorm")
fig = plt.figure("gH2AXnorm")

m = max(gH2AX)
for k in range(len(gH2AX)):
    gH2AXnorm[k] = 100 * gH2AX[k]/m


ax1 = fig.add_subplot(111)

ax1.scatter(t, gH2AXnorm, s=1, c='b', marker="s", label='gH2AXnorm')
xD = [ 23.90410959, 10,           5,           1.47260274, 0.71917808, 0.239726027, 0.136986301, 0.102739726 ]
yD = [ 10.55555556, 15.13888889, 26.52777778, 80,        100,          60,         20,           8.333333333 ]
ax1.scatter(xD, yD, s=20, c='r', marker="o", label='gH2AXnorm_dous')
xA = [  0.5308219178082174, 2.0376712328767113, 4.006849315068493, 11.952054794520548, 23.869863013698634 ]
yA = [ 99.99999999999997,  55.13888888888887,  30.27777777777777,  12.777777777777771,  1.6666666666666572 ]
ax1.scatter(xA, yA, s=30, c='g', marker="+", label='gH2AXnorm_Asaithamby')
plt.legend(loc='upper right');    

#plt.plot(t, gH2AXnorm)
plt.savefig("g-H2AX_norm.pdf")

plt.show()

with open(outputFile, 'w') as ofile:
    for m in range(len(gH2AX)):
        ofile.write(str(t[m])+ "\t" + str(gH2AXnorm[m])+"\n") 

#plt.xlim([1e1, 1e4])
#plt.ylim([1e-15, 3e-11])

print("\nCheck your folder for new images (pdf). Thank you for using molecularDNArepair.!\n")
