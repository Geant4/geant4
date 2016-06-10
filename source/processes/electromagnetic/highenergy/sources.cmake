#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emhighenergy
# Package: Geant4.src.G4processes.G4electromagnetic.G4emhighenergy
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66816 2013-01-12 16:13:54Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/muons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emhighenergy
    HEADERS
        G4AnnihiToMuPair.hh
        G4BetheBlochNoDeltaModel.hh
        G4BraggNoDeltaModel.hh
        G4GammaConversionToMuons.hh
        G4ICRU73NoDeltaModel.hh
        G4Vee2hadrons.hh
        G4ee2KChargedModel.hh
        G4ee2KNeutralModel.hh
        G4eeCrossSections.hh
        G4eeTo3PiModel.hh
        G4eeToHadrons.hh
        G4eeToHadronsModel.hh
        G4eeToHadronsMultiModel.hh
        G4eeToPGammaModel.hh
        G4eeToTwoPiModel.hh
        G4hBremsstrahlung.hh
        G4hBremsstrahlungModel.hh
        G4hPairProduction.hh
        G4hPairProductionModel.hh
        G4hhIonisation.hh
        G4mplIonisation.hh
        G4mplIonisationModel.hh
        G4mplIonisationWithDeltaModel.hh
    SOURCES
        G4AnnihiToMuPair.cc
        G4BetheBlochNoDeltaModel.cc
        G4BraggNoDeltaModel.cc
        G4GammaConversionToMuons.cc
        G4ICRU73NoDeltaModel.cc
        G4ee2KChargedModel.cc
        G4ee2KNeutralModel.cc
        G4eeCrossSections.cc
        G4eeTo3PiModel.cc
        G4eeToHadrons.cc
        G4eeToHadronsModel.cc
        G4eeToHadronsMultiModel.cc
        G4eeToPGammaModel.cc
        G4eeToTwoPiModel.cc
        G4hBremsstrahlung.cc
        G4hBremsstrahlungModel.cc
        G4hPairProduction.cc
        G4hPairProductionModel.cc
        G4hhIonisation.cc
        G4mplIonisation.cc
        G4mplIonisationModel.cc
        G4mplIonisationWithDeltaModel.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4emstandard
        G4emutils
        G4geometrymng
        G4globman
        G4hepnumerics
        G4intercoms
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4muons
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

