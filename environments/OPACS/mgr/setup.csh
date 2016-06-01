#  Installer, edit this file to set 
# correctly the product pathes : OPACS, 
# Geant4, CLHEP, RW, etc...
#
#  User, to use Geant4/environments/OPACS,
# execute this script with :
#   csh> source setup.csh
# 
#-----------------------------
#
# Execute OPACS setup script : 
# --------------------------
source /projects/OPACS/setup.csh
#
# Geant4 global path :
# ------------------
#setenv G4INSTALL /geant4/geant4beta
#
# If you have one, execute your Geant4 setup script :
# -------------------------------------------------
#source $G4INSTALL/mgr/setup.csh
#source $G4INSTALL/mgr/setsf.csh
#source $G4INSTALL/mgr/setiv.csh
#
# Template repository of Geant4 :
# -----------------------------
setenv OREPOSITORY $G4INSTALL/tmp/$G4SYSTEM/g4.ptrepository
#
# CLHEP path :
# ----------
setenv CLHEPROOT   $HOME/`uname`
#
# Rogue Wave :
# ----------
setenv RWROOT      $HOME/`uname`
#
# OPACS OPATH variable for finding, at run time, environments/OPACS palettes :
# -------------------------------------------------------------------------
setenv OPATH "$OPATH $G4INSTALL/environments/OPACS/usr"
#

