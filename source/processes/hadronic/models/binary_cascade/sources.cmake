#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_binary
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_binary
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_binary
    HEADERS
        G4Absorber.hh
        G4AntiProtonField.hh
        G4BCLateParticle.hh
        G4BinaryCascade.hh
        G4BinaryLightIonReaction.hh
        G4FieldPropagation.hh
        G4GeneratorPrecompoundInterface.hh
        G4KM_DummyField.hh
        G4KM_NucleonEqRhs.hh
        G4KM_OpticalEqRhs.hh
        G4KaonMinusField.hh
        G4KaonPlusField.hh
        G4KaonZeroField.hh
        G4NeutronField.hh
        G4PionMinusField.hh
        G4PionPlusField.hh
        G4PionZeroField.hh
        G4ProtonField.hh
        G4RKFieldIntegrator.hh
        G4RKPropagation.hh
        G4SigmaMinusField.hh
        G4SigmaPlusField.hh
        G4SigmaZeroField.hh
        G4VFieldPropagation.hh
        G4VKM_NuclearDensity.hh
        G4VNuclearField.hh
    SOURCES
        G4Absorber.cc
        G4AntiProtonField.cc
        G4BinaryCascade.cc
        G4BinaryLightIonReaction.cc
        G4FieldPropagation.cc
        G4GeneratorPrecompoundInterface.cc
        G4KM_NucleonEqRhs.cc
        G4KM_OpticalEqRhs.cc
        G4KaonMinusField.cc
        G4KaonPlusField.cc
        G4KaonZeroField.cc
        G4NeutronField.cc
        G4PionMinusField.cc
        G4PionPlusField.cc
        G4PionZeroField.cc
        G4ProtonField.cc
        G4RKFieldIntegrator.cc
        G4RKPropagation.cc
        G4SigmaMinusField.cc
        G4SigmaPlusField.cc
        G4SigmaZeroField.cc
        G4VFieldPropagation.cc
        G4VNuclearField.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_im_r_matrix
        G4had_mod_man
        G4had_mod_util
        G4had_preequ_exciton
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
        G4magneticfield
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

