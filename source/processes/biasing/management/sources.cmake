#------------------------------------------------------------------------------
# sources.cmake
# Module : G4biasing_mgt
# Package: Geant4.src.G4processes.G4biasing_mgt
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 94224 2015-11-09 08:30:58Z gcosmo $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/transportation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4biasing_mgt
    HEADERS
        G4VProcessPlacer.hh
        G4ProcessPlacer.hh
        G4BiasingAppliedCase.hh
        G4BiasingOperationManager.hh
        G4VBiasingInteractionLaw.hh
        G4VBiasingOperation.hh
        G4VBiasingOperator.hh
    SOURCES
        G4VProcessPlacer.cc
        G4ProcessPlacer.cc
        G4BiasingOperationManager.cc
        G4VBiasingOperation.cc
        G4VBiasingOperator.cc
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

