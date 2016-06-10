#------------------------------------------------------------------------------
# sources.cmake
# Module : G4leptons
# Package: Geant4.src.G4particles.G4leptons
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66892 2013-01-17 10:57:59Z gunter $
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
GEANT4_DEFINE_MODULE(NAME G4leptons
    HEADERS
        G4AntiNeutrinoE.hh
        G4AntiNeutrinoMu.hh
        G4AntiNeutrinoTau.hh
        G4Electron.hh
        G4LeptonConstructor.hh
        G4MuonMinus.hh
        G4MuonPlus.hh
        G4NeutrinoE.hh
        G4NeutrinoMu.hh
        G4NeutrinoTau.hh
        G4Positron.hh
        G4TauMinus.hh
        G4TauPlus.hh
    SOURCES
        G4AntiNeutrinoE.cc
        G4AntiNeutrinoMu.cc
        G4AntiNeutrinoTau.cc
        G4Electron.cc
        G4LeptonConstructor.cc
        G4MuonMinus.cc
        G4MuonPlus.cc
        G4NeutrinoE.cc
        G4NeutrinoMu.cc
        G4NeutrinoTau.cc
        G4Positron.cc
        G4TauMinus.cc
        G4TauPlus.cc
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

