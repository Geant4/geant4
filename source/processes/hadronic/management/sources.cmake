#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_mgt
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_mgt
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_mgt
    HEADERS
        G4EnergyRangeManager.hh
        G4HadLeadBias.hh
        G4HadronInelasticProcess.hh
        G4HadronicEPTestMessenger.hh
        G4HadronicProcess.hh
        G4HadronicProcessStore.hh
        G4HadronicProcessType.hh
        G4NoModelFound.hh
        G4VLeadingParticleBiasing.hh
    SOURCES
        G4EnergyRangeManager.cc
        G4HadLeadBias.cc
        G4HadronInelasticProcess.cc
        G4HadronicEPTestMessenger.cc
        G4HadronicProcess.cc
        G4HadronicProcessStore.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4hadronic_util
        G4hadronic_xsect
        G4intercoms
        G4ions
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
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

