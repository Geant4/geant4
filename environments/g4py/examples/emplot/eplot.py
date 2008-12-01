#!/usr/bin/python
# ==================================================================
# An example of ploting by EmCalculator
#
# Plotting photon cross sections and stopping power
# ==================================================================
from Geant4 import *
import g4py.ExN03pl
import g4py.emcalculator
import EmPlot

# initialize
EmPlot.Configure()

# user physics list
g4py.ExN03pl.Construct()

# target material
material= "G4_Cu"
EmPlot.SetMaterial(material)

# initialize G4 kernel
gRunManager.Initialize()
gRunManagerKernel.RunInitialization()

# energy
elist= []
for n in range(-3, 3):
  for i in range(10,99):
    elist.append(i/10.*10.**n *MeV)

# calculate stopping power
pname= "e-"
dedx_list= g4py.emcalculator.CalculateDEDX(pname, material, elist, 1)
xlist_tot=[]
xlist_ioni=[]
xlist_brems=[]

for x in dedx_list:
  xlist_tot.append((x[0], x[1]["tot"]/(MeV*cm2/g)))
  xlist_ioni.append((x[0], x[1]["ioni"]/(MeV*cm2/g)))
  xlist_brems.append((x[0], x[1]["brems"]/(MeV*cm2/g)))

# make plot
myCanvas= EmPlot.init_root()
aplot= EmPlot.make_plot(xlist_tot, pname+" Stopping Power ("+material+")",
                        "dE/dX (MeV cm^{2}/g)")
bplot= EmPlot.make_plot(xlist_ioni, "Stopping Power ("+material+")",
                        "dE/dX (MeV cm^{2}/g)", 1)
cplot= EmPlot.make_plot(xlist_brems, "Stopping Power ("+material+")",
                        "dE/dX (MeV cm^{2}/g)", 3)

