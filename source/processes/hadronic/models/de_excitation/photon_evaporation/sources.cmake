#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_photon_evaporation
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_photon_evaporation
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_photon_evaporation
    HEADERS
        G4GammaTransition.hh
        G4NeutronRadCapture.hh
        G4PhotonEvaporation.hh
        G4PolarizationTransition.hh
        G4VGammaTransition.hh
    SOURCES
        G4GammaTransition.cc
        G4NeutronRadCapture.cc
        G4PhotonEvaporation.cc
        G4PolarizationTransition.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4emlowenergy
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_deex_management
        G4hadronic_deex_util
        G4hadronic_util
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
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

