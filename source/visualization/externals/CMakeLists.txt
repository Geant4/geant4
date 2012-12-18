#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : 
# Package: Geant4.src.G4visualization.
#
# Intermediate level CMakeLists.txt - just process subdirectories
#
# Generated on : 24/9/2010
#
# $Id$
#
#------------------------------------------------------------------------------

#
# GL2PS built for GL and Inventor drivers
#
if(GEANT4_USE_OPENGL OR GEANT4_USE_INVENTOR)
    add_subdirectory(gl2ps)
endif()

#
# Geant4 internal ZLIB build optional
#
if(NOT GEANT4_USE_SYSTEM_ZLIB)
    add_subdirectory(zlib)
endif()

