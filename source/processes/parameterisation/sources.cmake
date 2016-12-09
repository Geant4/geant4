#------------------------------------------------------------------------------
# sources.cmake
# Module : G4parameterisation
# Package: Geant4.src.G4processes.G4parameterisation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 101152 2016-11-08 08:07:39Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4parameterisation
    HEADERS
        G4FastSimulationHelper.hh
        G4FastSimulationManager.hh
        G4FastSimulationManagerProcess.hh
        G4FastSimulationMessenger.hh
        G4FastSimulationProcessType.hh
        G4FastSimulationVector.hh
        G4FastSimulationVector.icc
        G4FastStep.hh
        G4FastStep.icc
        G4FastTrack.hh
        G4GlobalFastSimulationManager.hh
        G4VFastSimulationModel.hh
    SOURCES
        G4FastSimulationHelper.cc
        G4FastSimulationManager.cc
        G4FastSimulationManagerProcess.cc
        G4FastSimulationMessenger.cc
        G4FastStep.cc
        G4FastTrack.cc
        G4GlobalFastSimulationManager.cc
        G4VFastSimulationModel.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
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

