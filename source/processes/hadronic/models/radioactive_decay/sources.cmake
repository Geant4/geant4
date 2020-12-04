#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_radioactivedecay
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_radioactivedecay
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106689 2017-10-19 16:20:21Z dwright $
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4hadronic_radioactivedecay
    HEADERS
        G4AlphaDecay.hh
        G4BatemanParameters.hh
        G4BetaDecayCorrections.hh
        G4BetaDecayType.hh
        G4BetaMinusDecay.hh
        G4BetaPlusDecay.hh
        G4ECDecay.hh
        G4ITDecay.hh
        G4SFDecay.hh
        G4NeutronDecay.hh
        G4NuclearDecay.hh
        G4NucleusLimits.hh
        G4ProtonDecay.hh
        G4RadioactiveDecay.hh
        G4RadioactiveDecayBase.hh
        G4RadioactiveDecayBaseMessenger.hh
        G4Radioactivation.hh
        G4RadioactivationMessenger.hh
        G4RadioactiveDecayMode.hh
        G4RadioactiveDecayRatesToDaughter.hh
        G4RadioactiveDecayChainsFromParent.hh
        G4RadioactiveDecaymessenger.hh
        G4RadioactivityTable.hh
        G4TritonDecay.hh
        G4UIcmdWithNucleusLimits.hh
        G4UserLimitsForRD.hh
    SOURCES
        G4AlphaDecay.cc
        G4BatemanParameters.cc
        G4BetaDecayCorrections.cc
        G4BetaDecayType.cc
        G4BetaMinusDecay.cc
        G4BetaPlusDecay.cc
        G4ECDecay.cc
        G4ITDecay.cc
        G4SFDecay.cc
        G4NeutronDecay.cc
        G4NuclearDecay.cc
        G4NucleusLimits.cc
        G4ProtonDecay.cc
        G4RadioactiveDecay.cc
        G4RadioactiveDecayBase.cc
        G4RadioactiveDecayBaseMessenger.cc
        G4Radioactivation.cc
        G4RadioactivationMessenger.cc
        G4RadioactiveDecayMode.cc
        G4RadioactiveDecayRatesToDaughter.cc
        G4RadioactiveDecayChainsFromParent.cc
        G4RadioactiveDecaymessenger.cc
        G4RadioactivityTable.cc
        G4TritonDecay.cc
        G4UIcmdWithNucleusLimits.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4emutils
        G4geometrymng
        G4globman
        G4had_mod_util
        G4hadronic_deex_management
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4had_fission
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_util
        G4hadronic_xsect
        G4intercoms
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
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here
