#------------------------------------------------------------------------------
# sources.cmake
# Module : G4phys_ctor_em
# Package: Geant4.src.G4physicslists.G4physlist_ctors.G4physlist_ctor_em
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
GEANT4_DEFINE_MODULE(NAME G4phys_ctor_em
  HEADERS
    G4EmBuilder.hh
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
    G4EmLEPTSPhysics.hh
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
  SOURCES
    G4EmBuilder.cc
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
    G4EmLEPTSPhysics.cc
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
  GRANULAR_DEPENDENCIES
    G4baryons
    G4bosons
    G4cuts
    G4decay
    G4emdna-utils
    G4emdna-processes
    G4emdna-molman
    G4emdna-moltypes
    G4emdna-models
    G4emhighenergy
    G4emlowenergy
    G4emstandard
    G4emutils
    G4geometrymng
    G4globman
    G4had_mod_man
    G4hadronic_mgt
    G4hadronic_proc
    G4hadronic_xsect
    G4hadronic_util
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
    G4physlist_util
    G4procman
    G4run
    G4shortlived
    G4track
    G4transportation
    G4volumes
    G4xrays
  GLOBAL_DEPENDENCIES
    G4geometry
    G4global
    G4intercoms
    G4materials
    G4particles
    G4processes
    G4run
    G4track
  LINK_LIBRARIES
)

# List any source specific properties here

