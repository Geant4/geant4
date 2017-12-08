#------------------------------------------------------------------------------
# sources.cmake
# Module : G4biasing_imp
# Package: Geant4.src.G4processes.G4biasing_imp
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 105924 2017-08-29 12:28:09Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/biasing/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/biasing/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/transportation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4biasing_imp
    HEADERS
        G4GeometrySampler.hh
        G4ImportanceConfigurator.hh
        G4ImportanceProcess.hh
        G4PlaceOfAction.hh
        G4SamplingPostStepAction.hh
        G4VSampler.hh
        G4VSamplerConfigurator.hh
        G4WeightCutOffConfigurator.hh
        G4WeightCutOffProcess.hh
        G4WeightWindowConfigurator.hh
        G4WeightWindowProcess.hh
    SOURCES
        G4GeometrySampler.cc
        G4ImportanceConfigurator.cc
        G4ImportanceProcess.cc
        G4SamplingPostStepAction.cc
        G4VSampler.cc
        G4VSamplerConfigurator.cc
        G4WeightCutOffConfigurator.cc
        G4WeightCutOffProcess.cc
        G4WeightWindowConfigurator.cc
        G4WeightWindowProcess.cc
    GRANULAR_DEPENDENCIES
        G4cuts
        G4detector
        G4emutils
        G4geombias
        G4geometrymng
        G4globman
        G4hits
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4transportation
        G4volumes
	G4biasing_mgt
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

