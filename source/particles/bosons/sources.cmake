#------------------------------------------------------------------------------
# Module : G4bosons
# Package: Geant4.src.G4particles.G4bosons
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4bosons
  HEADERS
    G4BosonConstructor.hh
    G4ChargedGeantino.hh
    G4Gamma.hh
    G4Geantino.hh
    G4OpticalPhoton.hh
	G4PhononLong.hh
	G4PhononTransFast.hh
	G4PhononTransSlow.hh
    G4UnknownParticle.hh
  SOURCES
    G4BosonConstructor.cc
    G4ChargedGeantino.cc
    G4Gamma.cc
    G4Geantino.cc
    G4OpticalPhoton.cc
	G4PhononLong.cc
	G4PhononTransFast.cc
	G4PhononTransSlow.cc
    G4UnknownParticle.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4materials
    G4partman
  GLOBAL_DEPENDENCIES
    G4global
    G4materials
)

# List any source specific properties here
