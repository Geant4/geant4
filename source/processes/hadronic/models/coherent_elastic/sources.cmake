# - G4hadronic_coherent_elastic module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_coherent_elastic
  PUBLIC_HEADERS
    G4AntiNuclElastic.hh
    G4ChipsElasticModel.hh
    G4ChargeExchange.hh
    G4ChargeExchangeProcess.hh
    G4DiffuseElastic.hh
    G4DiffuseElasticV2.hh
    G4ElasticHadrNucleusHE.hh
    G4HadronElastic.hh
    G4LEHadronProtonElastic.hh
    G4LowEHadronElastic.hh
    G4hhElastic.hh
    G4LEnp.hh
    G4LEnpData.hh
    G4LEpp.hh
    G4LEppData.hh
    G4LMsdGenerator.hh
    G4NeutrinoElectronNcModel.hh
    G4NeutronElectronElModel.hh
    G4NuclNuclDiffuseElastic.hh
  SOURCES
    G4AntiNuclElastic.cc
    G4ChipsElasticModel.cc
    G4ChargeExchange.cc
    G4ChargeExchangeProcess.cc
    G4DiffuseElastic.cc
    G4DiffuseElasticV2.cc
    G4ElasticHadrNucleusHE.cc
    G4HadronElastic.cc
    G4LEHadronProtonElastic.cc
    G4LowEHadronElastic.cc
    G4hhElastic.cc
    G4LEnp.cc
    G4LEpp.cc
    G4LMsdGenerator.cc
    G4NeutrinoElectronNcModel.cc
    G4NeutronElectronElModel.cc
    G4NuclNuclDiffuseElastic.cc)

geant4_module_link_libraries(G4hadronic_coherent_elastic
  PUBLIC
    G4globman
    G4hadronic_mgt
    G4hadronic_util
    G4hadronic_xsect
    G4hepnumerics
    G4partman
    G4track
  PRIVATE
    G4baryons
    G4cuts
    G4hadronic_deex_handler
    G4heprandom
    G4ions
    G4leptons
    G4materials
    G4mesons)
