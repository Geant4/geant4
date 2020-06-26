#------------------------------------------------------------------------------
# Module : 
# Package: Geant4.src.G4global.
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4heprandom 
  HEADERS
    G4Poisson.hh
    G4QuickRand.hh
    G4RandomDirection.hh
    G4RandomTools.hh
	G4UniformRandPool.hh
    Randomize.hh
  SOURCES
    G4Poisson.cc
    G4UniformRandPool.cc
  GRANULAR_DEPENDENCIES
    G4globman
)

# List any source specific properties here
