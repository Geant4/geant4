#------------------------------------------------------------------------------
# Module : G4partadj
# Package: Geant4.src.G4particles.G4partadj
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4partadj
  HEADERS
    G4AdjointAlpha.hh
    G4AdjointDeuteron.hh
    G4AdjointElectron.hh
    G4AdjointElectronFI.hh
    G4AdjointGamma.hh
    G4AdjointGenericIon.hh
    G4AdjointHe3.hh
    G4AdjointIons.hh
    G4AdjointPositron.hh
    G4AdjointProton.hh
    G4AdjointTriton.hh
  SOURCES
    G4AdjointAlpha.cc
    G4AdjointDeuteron.cc
    G4AdjointElectron.cc
    G4AdjointElectronFI.cc
    G4AdjointGamma.cc
    G4AdjointGenericIon.cc
    G4AdjointHe3.cc
    G4AdjointIons.cc
    G4AdjointPositron.cc
    G4AdjointProton.cc
    G4AdjointTriton.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4materials
    G4partman
  GLOBAL_DEPENDENCIES
    G4global
    G4materials
)

# List any source specific properties here
