#------------------------------------------------------------------------------
# sources.cmake
# Module : G4muons
# Package: Geant4.src.G4processes.G4electromagnetic.G4muons
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 96156 2016-03-21 08:10:21Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4muons
    HEADERS
        G4EnergyLossForExtrapolator.hh
        G4ErrorEnergyLoss.hh
        G4MuBetheBlochModel.hh
        G4MuBremsstrahlung.hh
        G4MuBremsstrahlungModel.hh
        G4MuIonisation.hh
        G4MuMultipleScattering.hh
        G4MuPairProduction.hh
        G4MuPairProductionModel.hh
        G4TablesForExtrapolator.hh
        G4ePairProduction.hh
    SOURCES
        G4EnergyLossForExtrapolator.cc
        G4ErrorEnergyLoss.cc
        G4MuBetheBlochModel.cc
        G4MuBremsstrahlung.cc
        G4MuBremsstrahlungModel.cc
        G4MuIonisation.cc
        G4MuMultipleScattering.cc
        G4MuPairProduction.cc
        G4MuPairProductionModel.cc
        G4TablesForExtrapolator.cc
        G4ePairProduction.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4emstandard
        G4emutils
        G4geometrymng
        G4globman
        G4leptons
        G4materials
        G4mesons
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

