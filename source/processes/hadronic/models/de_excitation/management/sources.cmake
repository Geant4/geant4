#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_management
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_management
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_management
    HEADERS
        G4DeexParametersMessenger.hh
        G4DeexPrecoParameters.hh
        G4LevelManager.hh
        G4LevelReader.hh
        G4NuclearLevelData.hh
        G4NucLevel.hh
        G4VEmissionProbability.hh
        G4VEvaporationChannel.hh
        G4VEvaporationFactory.hh
    SOURCES
        G4DeexParametersMessenger.cc
        G4DeexPrecoParameters.cc
        G4LevelManager.cc
        G4LevelReader.cc
        G4NuclearLevelData.cc
        G4NucLevel.cc
        G4VEmissionProbability.cc
        G4VEvaporationChannel.cc
        G4VEvaporationFactory.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_util
        G4hepnumerics
        G4intercoms
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4track
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

