#!/usr/bin/env python
# coding: utf-8

# To run, change the links and definitions in this script and 
# type in the terminal containing the file 
# $python3 molecularDNAsurvival.py

# Authors: 
# K. Chatzipapas (*), D. Sakata, H. N. Tran, M. Dordevic, S. Incerti et al.
# (*) contact: k.chatzipapas@yahoo.com

# Authors' note : This is a very early version of this calculation technique. 
#                 Needs to optimized for each different application. An updated 
#                 version will be released in the future.


# Define filenames
outputFile = "molecularDNAsurvival.txt"
iRootFile  = "../molecular-dna.root"
print("\nInput file:",iRootFile,"\n")

# Define cell parameters
r3 = 7100*1e-09 * 2500*1e-09 * 7100*1e-09 # a * b * c   / Chromosome size, as defined in the mac file, but in meters. If sphere, a=b=c
NBP = 6405886128 # bp / Length of the DNA chain in Mbp
mass = 997 * 4 * 3.141592 * r3 / 3     # waterDensity * 4/3 * pi * r3 in kg

# Integration parameters
t0 = 0.0   # Starting time point (dimensionless)
t1 = 336   # Final time point (dimensionless)
dt = 0.001 # Intergration time step (dimensionless)
numpoints = int( (t1 - t0) / dt )

# Survival model Parameters
# THESE PARAMETERS HAVE TO BE DEFINED FOR PROPER RESULTS/ EXPERIMENTAL PARAMETERS, 
# AS DESCRIBED IN (Stewart RD. Radiat Res. 2001 https://pubmed.ncbi.nlm.nih.gov/11554848/ )
maxDose = 9
DR1     = 720      #Gy/h SF-dose
bp      = 1        #Gbp
#DR2   = 600.      #Gy/h SF-time
#DR3   = 12000.    #Gy/h FAR

# usual value of gamma 0.25 (changes the end value of the curve)
gamma = 0.19  #other published values: 0.189328    
# Lamb1 the end of the curve Lamb2 is connected with Eta
Lamb1 = 0.10  #other published values: 0.0125959   # 0.671  h-1    or 2.77 h-1 (15 min repair half-time)
Lamb2 = 1     #other published values: 33062.9     # 0.0439 h-1 (15.8-h repair half-time)    or 0.0616 h-1 (11.3 h repair half-time)
# 
Beta1 = 0     #other published values: 0.0193207   # 0.00152 
Beta2 = 0
# Defines the curvature as well as the end value of the curve
Eta   = 0.00005  #other published values: #7.50595e-06 # 1.18e-04 h-1  or 0.00042 h-1
# name of the cell (used in the output)
cell  = "test"

# ONLY ADVANCED USERS AFTER THIS POINT
#####################################################################################################################
import ROOT as root
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import threading

# Definition of the survival model
def SF_system( T_, DR_, Y_, Sig1_, Sig2_, gamma_, Lamb1_, Lamb2_, Beta1_, Beta2_, Eta_ ):

    global T 
    global DR 
    global Y  
    global Sig1  
    global Sig2  
    global gamma 
    global Lamb1 
    global Lamb2
    global Beta1 
    global Beta2
    global Eta  
    
    T  = T_
    DR = DR_
    Y  = Y_
    Sig1  = Sig1_
    Sig2  = Sig2_
    gamma = gamma_
    Lamb1 = Lamb1_
    Lamb2 = Lamb2_
    Beta1 = Beta1_
    Beta2 = Beta2_
    Eta   = Eta_


    return T, DR, Y, Sig1, Sig2, gamma, Lamb1, Lamb2, Beta1, Beta2, Eta

def SF_function( t, x, p ):
    if(t<=T):
        #print("t<=T")
        dx[0] = DR*Y*Sig1 - Lamb1*x[0] - Eta*x[0]*(x[0]+x[1])
        dx[1] = DR*Y*Sig2 - Lamb2*x[1] - Eta*x[1]*(x[0]+x[1])
        dx[2] = Beta1*Lamb1*x[0] + Beta2*Lamb2*x[1]+gamma*Eta*(x[0]+x[1])*(x[0]+x[1])
        dx[3] = (1.-Beta1)*Lamb1*x[0] + (1.-Beta2)*Lamb2*x[1]+(1.-gamma)*Eta*(x[0]+x[1])*(x[0]+x[1])
    else:
        dx[0] =  - Lamb1*x[0] - Eta*x[0]*(x[0]+x[1])
        dx[1] =  - Lamb2*x[1] - Eta*x[1]*(x[0]+x[1])
        dx[2] = Beta1*Lamb1*x[0] + Beta2*Lamb2*x[1]+gamma*Eta*(x[0]+x[1])*(x[0]+x[1])
        dx[3] = (1.-Beta1)*Lamb1*x[0] + (1.-Beta2)*Lamb2*x[1]+(1.-gamma)*Eta*(x[0]+x[1])*(x[0]+x[1])
    return dx


DSBBPID = []
repDSB = 0
irrDSB = 0
acc_edep = 0
dose = 0
i = 0
eVtoJ = 1.60218e-19


# Analyze root files
f = root.TFile(iRootFile)

fTree = f.Get("tuples/damage")
for e in fTree:
    if (e.TypeClassification=="DSB" or e.TypeClassification=="DSB+" or e.TypeClassification=="DSB++"):
        DSBBPID.append((e.Event, e.BasePair))
        i += 1

gTree = f.Get("tuples/classification")
for e in gTree:
    repDSB += e.DSB
    irrDSB += e.DSBp + 2*e.DSBpp
    
ffTree = f.Get("tuples/chromosome_hits")
for ev in ffTree:
    acc_edep += (ev.e_chromosome_kev + ev.e_dna_kev) *1e3  # eV    


# Calculate the total number of DSBs
totalDSB = repDSB + irrDSB
# Calculate the absorbed dose
dose = acc_edep * eVtoJ / mass    # Dose in Gy

#print ("r3=", r3)
#print("Dose :", dose, "\nTotal DSB/Gy :", totalDSB/dose, "\nRepairable DSBs/Gy :", repDSB/dose, "\nIrrepairable DSBs/Gy :", irrDSB/dose)

NBP_of_cell = 4682000000   # This is a constant of the model, as defined by Stewart in the TLK (Stewart RD. Radiat Res. 2001 https://pubmed.ncbi.nlm.nih.gov/11554848/ )
NBPcell = NBP/NBP_of_cell
totalDSBperGyCell = totalDSB/NBPcell/dose
#totalDSBperGyCell = totalDSB/dose
print ("Dose :", dose, "\nTotal DSB/Gy/Cell :", totalDSBperGyCell)

kFactor = 1  # This factor is used to normalize irrepairable damage to lower values than 1.
S1 = repDSB/NBPcell/dose/kFactor
S2 = irrDSB/NBPcell/dose/kFactor
print("Repairable DSBs/Gy/cell :", S1, "\nIrrepairable DSBs/Gy/cell :", S2, "\n\n")


# In[ ]:


T  = 0
DR = 0
Y  = 0
kozi = False
sum = [0]*(maxDose+1)
outputfile = [""] * (maxDose+1)

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
#t = [(t1-t0) * float(i) / (numpoints - 1) for i in range(numpoints)]
t = np.linspace(t0, t1, numpoints)

dx = np.zeros(4)
x  = np.zeros(4)
dsbs = 0.0
Sig1 = S1
Sig2 = S2

filename =  "txt/molecularDNAsurvival_"+cell+".dat"
outputfile1 = open(filename,"w")
outputfile1.write("Dose \t Survival Fraction\n")

for j in range(maxDose+1): 
    #j=1
    sol = []
    time = []
    dose = j 
    expEndTime = dose/DR1 #h
    
    #////double StopTime  = (expEndTime+ n/min(Lamb1,Lamb2));
        
    arguments = SF_system( expEndTime, DR1,bp,Sig1,Sig2,gamma,Lamb1,Lamb2,Beta1,Beta2,Eta )
    p = [ T, DR, Y, Sig1, Sig2, gamma, Lamb1, Lamb2, Beta1, Beta2, Eta ]

    Stepper = "dopri5"  # dopri5,vode,dop853
    
    if (j==0):
        solver = ode(SF_function).set_integrator(Stepper, nsteps=1)
        print("#*#*#*#*#*# Ignore this warning #*#*#*#*#*#")
    else:
        solver = ode(SF_function).set_integrator(Stepper, nsteps=numpoints)
    solver.set_initial_value( x, t0 ).set_f_params(p)
    #print(solver.t,solver.y)
        
    # Repeatedly call the `integrate` method to advance the
    # solution to time t+dt, and save the solution in solver.
    while solver.t < t1: # and solver.successful():   
        solver.integrate(solver.t+dt) 
        
    sum[j] = solver.y[2]
    print("Survival probability of cell ("+cell+") after",j,"Gy is", np.exp(-sum[j]) )
    outputfile1.write(str(dose) +"\t"+ str(np.exp(-(solver.y[2]))) +"\n")

c = np.linspace(0, maxDose, num=maxDose+1)
sum=np.array(sum)
outputfile1.close()

print("\nCalculation concluded. Check the final image and close it to end the program.!")

plt.yscale('log')
plt.xlim([0, 8])
plt.ylim([1e-2, 1])
plt.plot(c,np.exp(-sum[:]))
plt.show()

print("\nThank you for using molecularDNAsurvival.!\n")

