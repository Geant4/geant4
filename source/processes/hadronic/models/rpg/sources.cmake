#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_RPG
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_RPG
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_RPG
    HEADERS
        G4RPGAntiKZeroInelastic.hh
        G4RPGAntiLambdaInelastic.hh
        G4RPGAntiNeutronInelastic.hh
        G4RPGAntiOmegaMinusInelastic.hh
        G4RPGAntiProtonInelastic.hh
        G4RPGAntiSigmaMinusInelastic.hh
        G4RPGAntiSigmaPlusInelastic.hh
        G4RPGAntiXiMinusInelastic.hh
        G4RPGAntiXiZeroInelastic.hh
        G4RPGFragmentation.hh
        G4RPGInelastic.hh
        G4RPGKLongInelastic.hh
        G4RPGKMinusInelastic.hh
        G4RPGKPlusInelastic.hh
        G4RPGKShortInelastic.hh
        G4RPGKZeroInelastic.hh
        G4RPGLambdaInelastic.hh
        G4RPGNeutronInelastic.hh
        G4RPGNucleonInelastic.hh
        G4RPGOmegaMinusInelastic.hh
        G4RPGPiMinusInelastic.hh
        G4RPGPiPlusInelastic.hh
        G4RPGPionInelastic.hh
        G4RPGPionSuppression.hh
        G4RPGProtonInelastic.hh
        G4RPGReaction.hh
        G4RPGSigmaMinusInelastic.hh
        G4RPGSigmaPlusInelastic.hh
        G4RPGStrangeProduction.hh
        G4RPGTwoBody.hh
        G4RPGTwoCluster.hh
        G4RPGXiMinusInelastic.hh
        G4RPGXiZeroInelastic.hh
    SOURCES
        G4RPGAntiKZeroInelastic.cc
        G4RPGAntiLambdaInelastic.cc
        G4RPGAntiNeutronInelastic.cc
        G4RPGAntiOmegaMinusInelastic.cc
        G4RPGAntiProtonInelastic.cc
        G4RPGAntiSigmaMinusInelastic.cc
        G4RPGAntiSigmaPlusInelastic.cc
        G4RPGAntiXiMinusInelastic.cc
        G4RPGAntiXiZeroInelastic.cc
        G4RPGFragmentation.cc
        G4RPGInelastic.cc
        G4RPGKMinusInelastic.cc
        G4RPGKPlusInelastic.cc
        G4RPGKZeroInelastic.cc
        G4RPGLambdaInelastic.cc
        G4RPGNeutronInelastic.cc
        G4RPGNucleonInelastic.cc
        G4RPGOmegaMinusInelastic.cc
        G4RPGPiMinusInelastic.cc
        G4RPGPiPlusInelastic.cc
        G4RPGPionInelastic.cc
        G4RPGPionSuppression.cc
        G4RPGProtonInelastic.cc
        G4RPGReaction.cc
        G4RPGSigmaMinusInelastic.cc
        G4RPGSigmaPlusInelastic.cc
        G4RPGStrangeProduction.cc
        G4RPGTwoBody.cc
        G4RPGTwoCluster.cc
        G4RPGXiMinusInelastic.cc
        G4RPGXiZeroInelastic.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
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
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

