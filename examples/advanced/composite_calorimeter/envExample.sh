#------------------------------------------------------------------
# Bash-shell script to be run before building/executing this example.
# It has two parts: one for the Aida setup (for the histograms,
# ntuples, and more in general for the data analysis); and one 
# more specific for this test-beam example.
# Please notice that before running this script, you have to define 
# the usual Geant4 variables (in particular, the variable
# G4ANALYSIS_USE must be set to 1).
# For Aida, the variable PLATF is set for Linux redhat 7.2 with
# gcc 2.95.2 : if you are using a different platform and/or
# compiler, please change PLATF accordingly!
#------------------------------------------------------------------
#
# --- Aida setup ---
export LHCXXTOP=/afs/cern.ch/sw/lhcxx/
export LHCXXVERS=4.0.5
export PLATF=redhat72/gcc-2.95.2/
 . $LHCXXTOP/share/LHCXX/$LHCXXVERS/install/sharedstart.sh
#
# --- Specific setup for this test-beam example ---
export CCAL_CONFPATH=./dataconf
export CCAL_SENSITIVECONF=g4testbeamhcal96.conf
export CCAL_GEOMETRYCONF=testbeamhcal96.conf
export CCAL_GLOBALPATH=./dataglobal
export CCAL_GEOMPATH=./datageom
export CCAL_VISPATH=./datavis
# ---




