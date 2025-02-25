# - G4phys_ctor_em module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_em
  PUBLIC_HEADERS
    G4ChargedUnknownPhysics.hh
    G4EmBuilder.hh
    G4EmDNABuilder.hh
    G4EmDNAChemistry.hh
    G4EmDNAChemistry_option1.hh
    G4EmDNAChemistry_option2.hh
    G4EmDNAChemistry_option3.hh
    G4EmDNAPhysics.hh
    G4EmDNAPhysics_option1.hh
    G4EmDNAPhysics_option2.hh
    G4EmDNAPhysics_option3.hh
    G4EmDNAPhysics_option4.hh
    G4EmDNAPhysics_option5.hh
    G4EmDNAPhysics_option6.hh
    G4EmDNAPhysics_option7.hh
    G4EmDNAPhysics_option8.hh
    G4EmDNAPhysics_stationary.hh
    G4EmDNAPhysics_stationary_option2.hh
    G4EmDNAPhysics_stationary_option4.hh
    G4EmDNAPhysics_stationary_option6.hh
    G4EmDNAPhysicsActivator.hh
    G4EmLivermorePhysics.hh
    G4EmLivermorePolarizedPhysics.hh
    G4EmLowEPPhysics.hh
    G4EmModelActivator.hh
    G4EmParticleList.hh
    G4EmPenelopePhysics.hh
    G4EmStandardPhysics.hh
    G4EmStandardPhysicsGS.hh
    G4EmStandardPhysicsSS.hh
    G4EmStandardPhysicsWVI.hh
    G4EmStandardPhysics_option1.hh
    G4EmStandardPhysics_option2.hh
    G4EmStandardPhysics_option3.hh
    G4EmStandardPhysics_option4.hh
    G4GammaGeneralProcess.hh
    G4OpticalPhysics.hh
    G4ChemDissociationChannels.hh
    G4ChemDissociationChannels_option1.hh
  SOURCES
    G4ChargedUnknownPhysics.cc
    G4EmBuilder.cc
    G4EmDNABuilder.cc
    G4EmDNAChemistry.cc
    G4EmDNAChemistry_option1.cc
    G4EmDNAChemistry_option2.cc
    G4EmDNAChemistry_option3.cc
    G4EmDNAPhysics.cc
    G4EmDNAPhysics_option1.cc
    G4EmDNAPhysics_option2.cc
    G4EmDNAPhysics_option3.cc
    G4EmDNAPhysics_option4.cc
    G4EmDNAPhysics_option5.cc
    G4EmDNAPhysics_option6.cc
    G4EmDNAPhysics_option7.cc
    G4EmDNAPhysics_option8.cc
    G4EmDNAPhysics_stationary.cc
    G4EmDNAPhysics_stationary_option2.cc
    G4EmDNAPhysics_stationary_option4.cc
    G4EmDNAPhysics_stationary_option6.cc
    G4EmDNAPhysicsActivator.cc
    G4EmLivermorePhysics.cc
    G4EmLivermorePolarizedPhysics.cc
    G4EmLowEPPhysics.cc
    G4EmModelActivator.cc
    G4EmParticleList.cc
    G4EmPenelopePhysics.cc
    G4EmStandardPhysics.cc
    G4EmStandardPhysicsGS.cc
    G4EmStandardPhysicsSS.cc
    G4EmStandardPhysicsWVI.cc
    G4EmStandardPhysics_option1.cc
    G4EmStandardPhysics_option2.cc
    G4EmStandardPhysics_option3.cc
    G4EmStandardPhysics_option4.cc
    G4GammaGeneralProcess.cc
    G4OpticalPhysics.cc
    G4ChemDissociationChannels.cc
    G4ChemDissociationChannels_option1.cc)

geant4_module_link_libraries(G4phys_ctor_em
  PUBLIC
    G4emdna-utils
    G4emdna-processes
    G4emlowenergy
    G4emutils
    G4globman
    G4run
  PRIVATE
    G4baryons
    G4bosons
    G4cuts
    G4emdna-man
    G4emdna-models
    G4emdna-molman
    G4emdna-moltypes
    G4emhighenergy
    G4emstandard
    G4geometrymng
    G4hadronic_mgt
    G4hadronic_util
    G4ions
    G4leptons
    G4materials
    G4mesons
    G4muons
    G4optical
    G4partman
    G4phys_builders
    G4phys_ctor_factory
    G4physlist_util
    G4procman
    G4track
    G4transportation #NB Only for single enum in header
    G4xrays)
