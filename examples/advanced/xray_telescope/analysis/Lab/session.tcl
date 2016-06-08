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
delete storage
#
RioStorage storage g4osc.root READ
storage ls
#
delete hEnteringEnergy
delete hYZ
#
#  Get some histograms with the H1Get 
# builtin procedure defined in Lab/user/init.tcl :
#   H1DGet <variable> <storage> <SID>
H1DGet hEnteringEnergy storage "Entering energy"
#
# Plot the histo :
hEnteringEnergy vis
#
v set textContext "color black fontName TTF/couri"
v set histogramContext "color red modeling solid"
#
delete hYZ
H2DGet hYZ storage "YZ"
v set region 2
hYZ vis
