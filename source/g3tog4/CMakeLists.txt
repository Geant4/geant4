#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G3toG4
# Package: Geant4.src.G3toG4
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

# We also need to build the rztog4 application
# It's not a Mac Bundle or WIN32 app (i.e. no GUI!!)
# Will need cernlib libraries, but no FindPackage for that yet...
# Hmm, means it'll probably be UNIX only...
# NB: Doesn't build on Linux due to undefined references - so do
# need CERNLIB
#
#include_directories(include)
#add_executable(rztog4 
#    src/g3routines.F  
#    src/jshape.F  
#    src/rztog4.F
#    src/g3tog4.F      
#    src/mztog4.F  
#    src/tog4.F
#)

