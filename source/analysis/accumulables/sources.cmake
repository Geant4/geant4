#------------------------------------------------------------------------------
# Module : G4accumulables
# Package: Geant4.src.G4analysis.G4accumulables
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4accumulables
  HEADERS
    G4MergeMode.hh
    G4AccumulableManager.hh
    G4AccumulableManager.icc
    G4Accumulable.hh
    G4Accumulable.icc
    G4VAccumulable.hh
    G4VAccumulable.icc
  SOURCES
    G4MergeMode.cc
    G4AccumulableManager.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
)

# List any source specific properties here
