# - G4bosons module build definition

# Define the Geant4 Module.
geant4_add_module(G4bosons
  PUBLIC_HEADERS
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
    G4UnknownParticle.cc)

geant4_module_link_libraries(G4bosons PUBLIC G4partman G4globman)
