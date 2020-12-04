#------------------------------------------------------------------------------
# sources.cmake
# Module : G4phys_lists
# Package: Geant4.src.G4physicslists.G4physlist_ctors.G4physlist_ctor_hinelastic
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 10/01/2013
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4phys_ctor_hinelastic
    HEADERS
        G4HadronInelasticQBBC.hh
        G4HadronPhysicsFTF_BIC.hh
        G4HadronPhysicsFTFP_BERT.hh
        G4HadronPhysicsFTFP_BERT_HP.hh
        G4HadronPhysicsFTFP_BERT_TRV.hh
        G4HadronPhysicsFTFP_BERT_ATL.hh
        G4HadronPhysicsFTFQGSP_BERT.hh
	G4HadronPhysicsNuBeam.hh
        G4HadronPhysicsQGS_BIC.hh
        G4HadronPhysicsQGSP_BERT.hh
        G4HadronPhysicsQGSP_BERT_HP.hh
        G4HadronPhysicsQGSP_BIC.hh
        G4HadronPhysicsQGSP_BIC_HP.hh
        G4HadronPhysicsQGSP_FTFP_BERT.hh
        G4HadronPhysicsINCLXX.hh
        G4HadronPhysicsShielding.hh
        G4HadronPhysicsShieldingLEND.hh
        G4VHadronPhysics.hh
        G4HadronPhysicsQGSP_BIC_AllHP.hh
    SOURCES
        G4HadronInelasticQBBC.cc
        G4HadronPhysicsFTF_BIC.cc
        G4HadronPhysicsFTFP_BERT.cc
        G4HadronPhysicsFTFP_BERT_HP.cc
        G4HadronPhysicsFTFP_BERT_TRV.cc
        G4HadronPhysicsFTFP_BERT_ATL.cc
        G4HadronPhysicsFTFQGSP_BERT.cc
	G4HadronPhysicsNuBeam.cc
        G4HadronPhysicsQGS_BIC.cc
        G4HadronPhysicsQGSP_BERT.cc
        G4HadronPhysicsQGSP_BERT_HP.cc
        G4HadronPhysicsQGSP_BIC.cc
        G4HadronPhysicsQGSP_BIC_HP.cc
        G4HadronPhysicsQGSP_FTFP_BERT.cc
        G4HadronPhysicsINCLXX.cc
        G4HadronPhysicsShielding.cc
        G4HadronPhysicsShieldingLEND.cc
        G4VHadronPhysics.cc
        G4HadronPhysicsQGSP_BIC_AllHP.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4decay
        G4digits
        G4emhighenergy
        G4emlowenergy
        G4emstandard
        G4emutils
        G4event
        G4geometrymng
        G4globman
        G4had_fission
        G4had_im_r_matrix
        G4had_lend
        G4had_mod_man
        G4had_mod_util
        G4had_lept_nuclear
        G4had_preequ_exciton
        G4had_string_diff
        G4had_string_frag
        G4had_string_man
        G4had_theo_max
        G4hadronic_bert_cascade
        G4hadronic_binary
        G4hadronic_coherent_elastic
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_fission
        G4hadronic_deex_gem_evaporation
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_inclxx_interface
        G4hadronic_inclxx_physics
        G4hadronic_inclxx_utils
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_qgstring
        G4hadronic_radioactivedecay
        G4hadronic_stop
        G4hadronic_util
        G4hadronic_xsect
        G4hits
        G4intercoms
        G4ions
        G4leptons
        G4magneticfield
        G4materials
        G4mesons
        G4muons
        G4navigation
        G4optical
        G4partman
        G4phys_builders
        G4phys_ctor_factory
        G4procman
        G4run
        G4shortlived
        G4track
        G4tracking
        G4transportation
        G4volumes
        G4xrays
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4processes
        G4run
        G4track
        G4tracking
    LINK_LIBRARIES
)

# List any source specific properties here

