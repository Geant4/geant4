#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_hetcpp_utils
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_hetcpp.G4hadronic_hetcpp_utils
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_abla_abla
    HEADERS
        G4AblaDataDefs.hh
        G4AblaRandom.hh
        G4AblaInterface.hh
        G4AblaDataFile.hh
        G4Abla.hh
        G4AblaVirtualData.hh

    SOURCES
        G4AblaVirtualData.cc
        G4Abla.cc
        G4AblaRandom.cc
        G4AblaDataFile.cc
        G4AblaInterface.cc

    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4hadronic_mgt
        G4hadronic_util
        G4hadronic_xsect
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4track
        G4volumes
        G4intercoms
        G4had_preequ_exciton
        G4had_mod_man
        G4had_mod_util
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util

    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
        G4intercoms

    SOURCES_EXCLUDE_FORMAT
        G4Abla.cc

LINK_LIBRARIES)# List any source specific properties here
