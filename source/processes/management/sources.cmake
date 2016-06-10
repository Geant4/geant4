#------------------------------------------------------------------------------
# sources.cmake
# Module : G4procman
# Package: Geant4.src.G4processes.G4procman
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 73928 2013-09-17 08:00:50Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4procman
    HEADERS
        G4ParticleTypes.hh
        G4ProcTblElement.hh
        G4ProcTblElement.icc
        G4ProcessAttribute.hh
        G4ProcessManager.hh
        G4ProcessManager.icc
        G4ProcessManagerMessenger.hh
        G4ProcessTable.hh
        G4ProcessTable.icc
        G4ProcessTableMessenger.hh
        G4ProcessType.hh
        G4ProcessVector.hh
        G4ProcessVector.icc
        G4VContinuousDiscreteProcess.hh
        G4VContinuousProcess.hh
        G4VDiscreteProcess.hh
        G4VProcess.hh
        G4VRestContinuousDiscreteProcess.hh
        G4VRestContinuousProcess.hh
        G4VRestDiscreteProcess.hh
        G4VRestProcess.hh
        G4WrapperProcess.hh
    SOURCES
        G4ProcTblElement.cc
        G4ProcessAttribute.cc
        G4ProcessManager.cc
        G4ProcessManagerMessenger.cc
        G4ProcessTable.cc
        G4ProcessTableMessenger.cc
        G4ProcessVector.cc
        G4VContinuousDiscreteProcess.cc
        G4VContinuousProcess.cc
        G4VDiscreteProcess.cc
        G4VProcess.cc
        G4VRestContinuousDiscreteProcess.cc
        G4VRestContinuousProcess.cc
        G4VRestDiscreteProcess.cc
        G4VRestProcess.cc
        G4WrapperProcess.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4intercoms
        G4materials
        G4partman
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

