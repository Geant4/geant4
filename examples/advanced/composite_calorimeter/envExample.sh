#------------------------------------------------------------------
# Bash-shell script to be run before building/executing this example.
#------------------------------------------------------------------
#
# --- Geant4 specific ---
export G4ANALYSIS_USE=1
#
# --- Aida / PI ---
eval `aida-config --runtime sh`
#
# --- Specific setup for this test-beam example ---
export CCAL_CONFPATH=./dataconf
export CCAL_SENSITIVECONF=g4testbeamhcal96.conf
export CCAL_GEOMETRYCONF=testbeamhcal96.conf
export CCAL_GLOBALPATH=./dataglobal
export CCAL_GEOMPATH=./datageom
export CCAL_VISPATH=./datavis
# ---




