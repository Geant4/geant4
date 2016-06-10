#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_multifragmentation
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_multifragmentation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88406 2015-02-18 09:13:29Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_multifragmentation
    HEADERS
        G4Solver.hh
        G4Solver.icc
        G4StatMF.hh
        G4StatMFChannel.hh
        G4StatMFFragment.hh
        G4StatMFMacroBiNucleon.hh
        G4StatMFMacroCanonical.hh
        G4StatMFMacroChemicalPotential.hh
        G4StatMFMacroMultiNucleon.hh
        G4StatMFMacroMultiplicity.hh
        G4StatMFMacroNucleon.hh
        G4StatMFMacroTemperature.hh
        G4StatMFMacroTetraNucleon.hh
        G4StatMFMacroTriNucleon.hh
        G4StatMFMicroCanonical.hh
        G4StatMFMicroManager.hh
        G4StatMFMicroPartition.hh
        G4StatMFParameters.hh
        G4VMultiFragmentation.hh
        G4VStatMFEnsemble.hh
        G4VStatMFMacroCluster.hh
    SOURCES
        G4Solver.cc
        G4StatMF.cc
        G4StatMFChannel.cc
        G4StatMFFragment.cc
        G4StatMFMacroBiNucleon.cc
        G4StatMFMacroCanonical.cc
        G4StatMFMacroChemicalPotential.cc
        G4StatMFMacroMultiNucleon.cc
        G4StatMFMacroMultiplicity.cc
        G4StatMFMacroNucleon.cc
        G4StatMFMacroTemperature.cc
        G4StatMFMacroTetraNucleon.cc
        G4StatMFMacroTriNucleon.cc
        G4StatMFMicroCanonical.cc
        G4StatMFMicroManager.cc
        G4StatMFMicroPartition.cc
        G4StatMFParameters.cc
        G4VMultiFragmentation.cc
        G4VStatMFEnsemble.cc
        G4VStatMFMacroCluster.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
        G4hadronic_mgt
        G4hadronic_util
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4track
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

