#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_string_man
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_string.G4had_string_man
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
GEANT4_DEFINE_MODULE(NAME G4had_string_man
    HEADERS
        G4InteractionContent.hh
        G4StringModel.hh
        G4VParticipants.hh
        G4VPartonStringModel.hh
        G4VSplitableHadron.hh
        G4VStringFragmentation.hh
    SOURCES
        G4InteractionContent.cc
        G4StringModel.cc
        G4VParticipants.cc
        G4VPartonStringModel.cc
        G4VSplitableHadron.cc
        G4VStringFragmentation.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_mgt
        G4hadronic_util
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4shortlived
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

