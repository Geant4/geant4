#------------------------------------------------------------------
# C-shell script to be run before building/executing this example.
# It has two parts: one for the Aida setup (for the histograms,
# ntuples, and more in general for the data analysis); and one 
# more specific for this test-beam example.
# Please notice that before running this script, you have to define 
# the usual Geant4 variables (in particular, the variable
# G4ANALYSIS_USE must be set to 1).
#------------------------------------------------------------------
#
### # --- Aida 2.2 , Anaphe 4.0.5 ---
### setenv LHCXXTOP   /afs/cern.ch/sw/lhcxx/
### setenv LHCXXVERS  4.0.5
### setenv PLATF      redhat72/gcc-2.95.2/
### source $LHCXXTOP/share/LHCXX/$LHCXXVERS/install/sharedstart.csh
#
# --- Aida 3.0 , Anaphe 5.0.1 ---
#source /afs/cern.ch/sw/lhcxx/share/LHCXX/5.0.1/scripts/setupAnaphe.csh
#setenv PATH ${PATH}:/afs/cern.ch/sw/lhcxx/share/LHCXX/5.0.1/scripts/
#
# --- Specific setup for this test-beam example ---
setenv CCAL_CONFPATH       ./dataconf
setenv CCAL_SENSITIVECONF  g4testbeamhcal96.conf
setenv CCAL_GEOMETRYCONF   testbeamhcal96.conf
setenv CCAL_GLOBALPATH     ./dataglobal
setenv CCAL_GEOMPATH       ./datageom
setenv CCAL_VISPATH        ./datavis
# ---




