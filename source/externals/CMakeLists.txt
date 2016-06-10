#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4externals
# Package: Geant4.src.G4externals
#
# Intermediate level CMakeLists.txt - just process subdirectories
#
# Created on : 02/06/2011
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# We only need internal CLHEP if it was requested
#
if(NOT GEANT4_USE_SYSTEM_CLHEP)
    add_subdirectory(clhep)
endif()

#----------------------------------------------------------------------------
# We only need internal expat if it was requested
#
if(NOT GEANT4_USE_SYSTEM_EXPAT)
    add_subdirectory(expat)
endif()

#----------------------------------------------------------------------------
# Geant4 internal ZLIB build optional
#
if(NOT GEANT4_USE_SYSTEM_ZLIB)
    add_subdirectory(zlib)
endif()
