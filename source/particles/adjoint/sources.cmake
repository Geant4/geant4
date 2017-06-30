#------------------------------------------------------------------------------
# sources.cmake
# Module : G4partadj
# Package: Geant4.src.G4particles.G4partadj
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 102170 2017-01-09 13:16:15Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4partadj
    HEADERS
        G4AdjointAlpha.hh
        G4AdjointDeuteron.hh
        G4AdjointElectron.hh
        G4AdjointElectronFI.hh
        G4AdjointGamma.hh
        G4AdjointGenericIon.hh
        G4AdjointHe3.hh
        G4AdjointIons.hh
        G4AdjointPositron.hh
        G4AdjointProton.hh
        G4AdjointTriton.hh
    SOURCES
        G4AdjointAlpha.cc
        G4AdjointDeuteron.cc
        G4AdjointElectron.cc
        G4AdjointElectronFI.cc
        G4AdjointGamma.cc
        G4AdjointGenericIon.cc
        G4AdjointHe3.cc
        G4AdjointIons.cc
        G4AdjointPositron.cc
        G4AdjointProton.cc
        G4AdjointTriton.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4materials
        G4partman
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
    LINK_LIBRARIES
)

# List any source specific properties here

