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
source /lal/OPACS/v3/setup.csh
#
# Geant4 global path :
# ------------------
#setenv G4INSTALL /geant4/dev/geant4
#
# If you have one, execute your Geant4 setup script :
# -------------------------------------------------
#  Mainly it should define G4INSTALL, G4SYSTEM, G4LIB, G4WORKDIR 
# needed by the mgr/Makefile, Config.mk make files :
#source $G4INSTALL/mgr/setup.csh
#source $G4INSTALL/mgr/setsf.csh
#source $G4INSTALL/mgr/setiv.csh
#
# CLHEP, Rogue Wave path :
# ----------------------
if ( `uname -n` == "aleph" ) then
setenv CLHEPROOT /lal/CLHEP/1.4/HP-UX-aCC
setenv RWLIBS    "-lrwtool"
endif
if ( `uname -n` == "asc" ) then
setenv CLHEPROOT /lal/CLHEP/1.4/OSF1-cxx
setenv RWROOT    /lal/rogue/6.1/OSF1-cxx
setenv RWFLAGS   "-I${RWROOT}"
setenv RWLIBS    "-L${RWROOT}/lib -lrwtool"
endif
if ( `uname -n` == "lx1.lal.in2p3.fr" ) then
setenv CLHEPROOT /lal/CLHEP/1.4/Linux-gxx
setenv RWROOT    /lal/rogue/6.1/Linux-gxx
setenv RWFLAGS   "-I${RWROOT}"
setenv RWLIBS    "-L${RWROOT}/lib -lrwtool"
endif
if ( `uname -n` == "papou1" ) then
setenv G4WORKDIR /geant4/dev
setenv G4INSTALL $G4WORKDIR/geant4
setenv G4LIB     $G4WORKDIR/lib
setenv G4SYSTEM  SunOS-CC
setenv CLHEPROOT /lal/CLHEP/1.4/SunOS-CC
setenv RWROOT    /lal/rogue/6.1/SunOS-CC
setenv RWFLAGS   "-I${RWROOT}"
setenv RWLIBS    "-L${RWROOT}/lib -lrwtool"
endif
#
# Directory where to find Geant4 templates (to be able to link) :
# -------------------------------------------------------------
setenv OREPOSITORY $G4WORKDIR/tmp/$G4SYSTEM/g4.ptrepository
#
# OPACS OPATH variable for finding, at run time, environments/OPACS palettes :
# -------------------------------------------------------------------------
setenv OPATH "$OPATH $G4INSTALL/environments/OPACS/usr"
#

