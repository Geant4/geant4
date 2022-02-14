# - G4phonon module build definition

# Define the Geant4 Module.
geant4_add_module(G4solidstate_phonon
  PUBLIC_HEADERS
	G4LatticeManager.hh
	G4LatticeReader.hh
	G4PhononDownconversion.hh
	G4PhononPolarization.hh
	G4PhononReflection.hh
	G4PhononScattering.hh
	G4PhononTrackMap.hh
	G4VPhononProcess.hh
  SOURCES
	G4LatticeManager.cc
	G4LatticeReader.cc
	G4PhononDownconversion.cc
	G4PhononPolarization.cc
	G4PhononReflection.cc
	G4PhononScattering.cc
	G4PhononTrackMap.cc
	G4VPhononProcess.cc)

geant4_module_link_libraries(G4solidstate_phonon
  PUBLIC
    G4globman
    G4procman
  PRIVATE
    G4bosons
    G4geometrymng
    G4heprandom
    G4materials
    G4partman
    G4track)
