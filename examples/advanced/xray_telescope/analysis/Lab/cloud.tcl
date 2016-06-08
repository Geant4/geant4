#
# ********************************************************************
# * DISCLAIMER                                                       *
# *                                                                  *
# * The following disclaimer summarizes all the specific disclaimers *
# * of contributors to this software. The specific disclaimers,which *
# * govern, are listed with their locations in:                      *
# *   http://cern.ch/geant4/license                                  *
# *                                                                  *
# * Neither the authors of this software system, nor their employing *
# * institutes,nor the agencies providing financial support for this *
# * work  make  any representation or  warranty, express or implied, *
# * regarding  this  software system or assume any liability for its *
# * use.                                                             *
# *                                                                  *
# * This  code  implementation is the  intellectual property  of the *
# * GEANT4 collaboration.                                            *
# * By copying,  distributing  or modifying the Program (or any work *
# * based  on  the Program)  you indicate  your  acceptance of  this *
# * statement, and all its terms.                                    *
# ********************************************************************
#
#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
delete v
#
# Get a SWIG handler over the current viewer :
IViewer v -this [ui getCurrentViewer]
v reset page
#
# 2x2 regions :
v set page "1 2"
v set pageTitle "Xray telescope"
#
#/////////////////////////////////////////////////////////////////////////////
#delete storage
#
#RioStorage storage g4osc.root READ
#storage ls
#
#delete hEnteringEnergy
#
#  Get some histograms with the H1Get 
# builtin procedure defined in Lab/user/init.tcl :
#   H1DGet <variable> <storage> <SID>
#H1DGet hEnteringEnergy storage "Entering energy"
#
# Plot the histo :
#hEnteringEnergy vis
#
#v set textContext "color black fontName TTF/couri"
#v set histogramContext "color red modeling solid"
#
#/////////////////////////////////////////////////////////////////////////////
delete s2
delete tuple
delete t_energy
#
v set region 1
#
RioStorage s2 XrayTel.root READ
RioTuple tuple [& s2] "XrayTel"
tuple print
# 
H1D t_energy "Entering energy" 100 0 0.002
set cuts ""
t_energy project [& tuple] EnteringEnergy $cuts
t_energy vis
unset cuts
#
# Fitting :
#delete fitGauss
#
#Gaussian fitGauss 100 0.0008 0.0012
#fitGauss fit [& t_energy]
#
#fitGauss vis
#
#/////////////////////////////////////////////////////////////////////////////
delete yz
#
v set region 2
#
Cloud2D yz YZ
set cuts ""
yz project [& tuple] y z $cuts
yz vis
unset cuts
#
v set cloudContext "modeling cloud color blue markerStyle plus markerSize 10"
#
