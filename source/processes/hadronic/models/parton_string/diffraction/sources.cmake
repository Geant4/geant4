#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_string_diff
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_string.G4had_string_diff
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
GEANT4_DEFINE_MODULE(NAME G4had_string_diff
    HEADERS
        G4DiffractiveExcitation.hh
        G4DiffractiveSplitableHadron.hh
        G4ElasticHNScattering.hh
        G4FTFAnnihilation.hh
        G4FTFModel.hh
        G4FTFParameters.hh
        G4FTFParticipants.hh
    SOURCES
        G4DiffractiveExcitation.cc
        G4DiffractiveSplitableHadron.cc
        G4ElasticHNScattering.cc
        G4FTFAnnihilation.cc
        G4FTFModel.cc
        G4FTFParameters.cc
        G4FTFParticipants.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4had_string_diff
        G4had_string_frag
        G4had_string_man
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

