#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////////
#
# Example of tcl script using a root file.
#
#/////////////////////////////////////////////////////////////////////////////
# If reexecuted :
delete storage
delete v
#
# Get a SWIG handler over the current viewer :
IViewer v -this [ui getCurrentViewer]
v reset page
#
# 2x2 regions :
#v set page "2 2"
v set pageTitle "G4 analysis with Open Scientist"
#
RioStorage storage g4osc.root READ
storage ls
#
delete EAbs
#
#storage cd histograms
#storage ls
# 
#  Get some histograms with the H1Get 
# builtin procedure defined in Lab/user/init.tcl :
#  In the below the first "EAbs" is a variable and
# the second is the name of the object in the storage
# (the SID, storage identifier) :
#   H1DGet <variable> <storage> <SID>
H1DGet EAbs storage EAbs
#
# Plot the histo :
EAbs vis
#
v set textContext "color black fontName TTF/couri"
v set histogramContext "color red modeling solid"
#
# Fit :
#delete fitExp
#
#Exponential fitExp 0. 1.
#fitExp fit [& EAbs]
#
#fitExp vis
#
