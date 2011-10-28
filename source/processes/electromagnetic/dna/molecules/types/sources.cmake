#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emlowenergy
# Package: Geant4.src.G4processes.G4electromagnetic.G4emlowenergy
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.4 2010-11-15 08:24:43 gcosmo Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/track/include)
include_directories(${CMAKE_SOURCE_DIR}/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/materials/include)
#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emdna-moltypes
    HEADERS
        G4Electron_aq.hh
        G4H2.hh
        G4H2O2.hh
        G4H2O.hh
        G4H3O.hh
        G4Hydrogen.hh
        G4OH.hh
    SOURCES
        G4Electron_aq.cc
        G4H2.cc
        G4H2O2.cc
        G4H2O.cc
        G4H3O.cc
        G4Hydrogen.cc
        G4OH.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4heprandom
        G4materials
        G4partman
        G4track
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

