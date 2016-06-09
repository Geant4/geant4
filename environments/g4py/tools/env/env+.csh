#$Id: env+.csh,v 1.2 2006-08-08 09:32:52 kmura Exp $
# ======================================================================
#  [csh version]
#  Additonal setup script for co-building global shared libraries
#
#  In dafault situation, Geant4 creates and uses "granular static libraries".
#  But, Python interface requires "global shared library".
#  This small setup script helps co-build global shared libraries,
#  in addition to an existing granular static environment.
#
#  After sourcing "env.csh", which is created by the "Configure" script,
#  source this file additionally, then "make global".
#  Global shared library staffs will be created in another directory.
# ======================================================================

if (! $?G4INSTALL) then
  echo "*** G4INSTALL is not defined! Please source env.csh first."
  exit
endif

# Recompilation is necessary because position-independent codes (-fPIC option) 
# shoud be required. The temporal directory should be different from 
# the existing one.
setenv G4TMP $G4INSTALL/tmp-slib

# Under this directory, new global shared library staffs will be created.
setenv G4LIB $G4INSTALL/slib

setenv G4LIB_BUILD_SHARED 1

