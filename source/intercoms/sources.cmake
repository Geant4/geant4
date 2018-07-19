#------------------------------------------------------------------------------
# sources.cmake
# Module : G4intercoms
# Package: Geant4.src.G4intercoms
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106172 2017-09-15 13:03:57Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4intercoms
    HEADERS
        G4AnyMethod.hh
        G4AnyType.hh
        G4GenericMessenger.hh
        G4LocalThreadCoutMessenger.hh
        G4UIaliasList.hh
        G4UIbatch.hh
        G4UIbridge.hh
        G4UIcmdWith3Vector.hh
        G4UIcmdWith3VectorAndUnit.hh
        G4UIcmdWithABool.hh
        G4UIcmdWithADouble.hh
        G4UIcmdWithADoubleAndUnit.hh
        G4UIcmdWithAString.hh
        G4UIcmdWithAnInteger.hh
        G4UIcmdWithoutParameter.hh
        G4UIcommand.hh
        G4UIcommandStatus.hh
        G4UIcommandTree.hh
        G4UIcontrolMessenger.hh
        G4UIdirectory.hh
        G4UImanager.hh
        G4UImessenger.hh
        G4UIparameter.hh
        G4UIsession.hh
        G4UItokenNum.hh
        G4UnitsMessenger.hh
        G4VFlavoredParallelWorld.hh
        G4VGlobalFastSimulationManager.hh
        icomsdefs.hh
    SOURCES
        G4LocalThreadCoutMessenger.cc
        G4GenericMessenger.cc
        G4UIaliasList.cc
        G4UIbatch.cc
        G4UIbridge.cc
        G4UIcmdWith3Vector.cc
        G4UIcmdWith3VectorAndUnit.cc
        G4UIcmdWithABool.cc
        G4UIcmdWithADouble.cc
        G4UIcmdWithADoubleAndUnit.cc
        G4UIcmdWithAString.cc
        G4UIcmdWithAnInteger.cc
        G4UIcmdWithoutParameter.cc
        G4UIcommand.cc
        G4UIcommandTree.cc
        G4UIcontrolMessenger.cc
        G4UIdirectory.cc
        G4UImanager.cc
        G4UImessenger.cc
        G4UIparameter.cc
        G4UIsession.cc
        G4UnitsMessenger.cc
        G4VGlobalFastSimulationManager.cc
    GRANULAR_DEPENDENCIES
        G4globman
    GLOBAL_DEPENDENCIES
        G4global
    LINK_LIBRARIES
)

# List any source specific properties here

