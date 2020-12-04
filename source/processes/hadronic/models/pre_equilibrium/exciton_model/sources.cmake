#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_preequ_exciton
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_preequ.G4had_preequ_exciton
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
GEANT4_DEFINE_MODULE(NAME G4had_preequ_exciton
    HEADERS
        G4GNASHTransitions.hh
        G4HETCAlpha.hh
        G4HETCChargedFragment.hh
        G4HETCDeuteron.hh
        G4HETCEmissionFactory.hh
        G4HETCFragment.hh
        G4HETCFragment.icc
        G4HETCHe3.hh
        G4HETCNeutron.hh
        G4HETCProton.hh
        G4HETCTriton.hh
        G4LowEGammaNuclearModel.hh
        G4LowEIonFragmentation.hh
        G4PreCompoundAlpha.hh
        G4PreCompoundDeuteron.hh
        G4PreCompoundEmission.hh
        G4PreCompoundEmissionFactory.hh
        G4PreCompoundFragment.hh
        G4PreCompoundFragmentVector.hh
        G4PreCompoundHe3.hh
        G4PreCompoundIon.hh
        G4PreCompoundModel.hh
        G4PreCompoundNeutron.hh
        G4PreCompoundNucleon.hh
        G4PreCompoundProton.hh
        G4PreCompoundTransitions.hh
        G4PreCompoundTriton.hh
        G4VPreCompoundEmissionFactory.hh
        G4VPreCompoundFragment.hh
        G4VPreCompoundFragment.icc
        G4VPreCompoundTransitions.hh
    SOURCES
        G4GNASHTransitions.cc
        G4HETCAlpha.cc
        G4HETCChargedFragment.cc
        G4HETCDeuteron.cc
        G4HETCEmissionFactory.cc
        G4HETCFragment.cc
        G4HETCHe3.cc
        G4HETCNeutron.cc
        G4HETCProton.cc
        G4HETCTriton.cc
        G4LowEGammaNuclearModel.cc
        G4LowEIonFragmentation.cc
        G4PreCompoundAlpha.cc
        G4PreCompoundDeuteron.cc
        G4PreCompoundEmission.cc
        G4PreCompoundEmissionFactory.cc
        G4PreCompoundFragment.cc
        G4PreCompoundFragmentVector.cc
        G4PreCompoundHe3.cc
        G4PreCompoundIon.cc
        G4PreCompoundModel.cc
        G4PreCompoundNeutron.cc
        G4PreCompoundNucleon.cc
        G4PreCompoundProton.cc
        G4PreCompoundTransitions.cc
        G4PreCompoundTriton.cc
        G4VPreCompoundEmissionFactory.cc
        G4VPreCompoundFragment.cc
        G4VPreCompoundTransitions.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_util
        G4hadronic_xsect
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

