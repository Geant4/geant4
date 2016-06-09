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
material= "G4_Pb"
EmPlot.SetMaterial(material)

# initialize G4 kernel
gRunManager.Initialize()
gRunManagerKernel.RunInitialization()

# energy
elist= []
for n in range(-3, 4):
  for i in range(10,99):
    elist.append(i/10.*10.**n *MeV)


# calculate cross sections
xsection_list= g4py.emcalculator.CalculatePhotonCrossSection(material, elist, 1)
xlist_tot=[]
xlist_comp=[]
xlist_pe=[]
xlist_conv=[]
for x in xsection_list:
  xlist_tot.append((x[0]/MeV, x[1]["tot"]/(cm2/g)))
  xlist_comp.append((x[0]/MeV, x[1]["compt"]/(cm2/g)))
  xlist_pe.append((x[0]/MeV, x[1]["phot"]/(cm2/g)))
  xlist_conv.append((x[0]/MeV, x[1]["conv"]/(cm2/g)))

# make a plot
myCanvas= EmPlot.init_root()
aplot= EmPlot.make_plot(xlist_tot,  "Photon Cross Section ("+material+")",
                        "Cross Section (cm^{2}/g)")
bplot= EmPlot.make_plot(xlist_comp, "Photon Cross Section ("+material+")",
                        "Cross Section (cm^{2}/g)", 1)
cplot= EmPlot.make_plot(xlist_pe,   "Photon Cross Section ("+material+")",
                        "Cross Section (cm^{2}/g)", 7)
dplot= EmPlot.make_plot(xlist_conv, "Photon Cross Section ("+material+")",
                        "Cross Section (cm^{2}/g)", 3)

