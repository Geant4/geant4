#------------------------------------------------------------------
# C-shell script to be run before building/executing this example.
#------------------------------------------------------------------
#
# --- Geant4 specific ---
setenv G4ANALYSIS_USE 1
#
# --- Aida / PI ---
eval `aida-config --runtime csh`
#
# --- Specific setup for this test-beam example ---
setenv CCAL_CONFPATH       ./dataconf
setenv CCAL_SENSITIVECONF  g4testbeamhcal96.conf
setenv CCAL_GEOMETRYCONF   testbeamhcal96.conf
setenv CCAL_GLOBALPATH     ./dataglobal
setenv CCAL_GEOMPATH       ./datageom
setenv CCAL_VISPATH        ./datavis
# ---




