#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_mod_man
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_mod_man
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
GEANT4_DEFINE_MODULE(NAME G4had_mod_man
    HEADERS
        G4HadronicInteraction.hh
        G4HadronicInteractionRegistry.hh
        G4V3DNucleus.hh
        G4VCrossSectionBase.hh
        G4VHighEnergyGenerator.hh
        G4VIntraNuclearTransportModel.hh
        G4VKineticNucleon.hh
        G4VNuclearDensity.hh
        G4VPreCompoundModel.hh
    SOURCES
        G4HadronicInteraction.cc
        G4HadronicInteractionRegistry.cc
        G4V3DNucleus.cc
        G4VCrossSectionBase.cc
        G4VHighEnergyGenerator.cc
        G4VIntraNuclearTransportModel.cc
        G4VKineticNucleon.cc
        G4VNuclearDensity.cc
        G4VPreCompoundModel.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4hadronic_util
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
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

