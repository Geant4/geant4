# - G4heprandom module build definition

geant4_add_module(G4heprandom
  PUBLIC_HEADERS
    G4Poisson.hh
    G4QuickRand.hh
    G4RandomDirection.hh
    G4RandomTools.hh
	  G4UniformRandPool.hh
    Randomize.hh
  SOURCES
    G4Poisson.cc
    G4UniformRandPool.cc)

geant4_module_link_libraries(G4heprandom PUBLIC G4globman)
