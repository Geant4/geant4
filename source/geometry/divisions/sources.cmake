#------------------------------------------------------------------------------
# Module : G4geomdivision
# Package: Geant4.src.G4geometry.G4geomdivision
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4geomdivision
  HEADERS
    G4PVDivision.hh
    G4PVDivisionFactory.hh
    G4ParameterisationBox.hh
    G4ParameterisationCons.hh
    G4ParameterisationPara.hh
    G4ParameterisationPolycone.hh
    G4ParameterisationPolyhedra.hh
    G4ParameterisationTrd.hh
    G4ParameterisationTubs.hh
    G4ReplicatedSlice.hh
    G4VDivisionParameterisation.hh
    G4VDivisionParameterisation.icc
  SOURCES
    G4PVDivision.cc
    G4PVDivisionFactory.cc
    G4ParameterisationBox.cc
    G4ParameterisationCons.cc
    G4ParameterisationPara.cc
    G4ParameterisationPolycone.cc
    G4ParameterisationPolyhedra.cc
    G4ParameterisationTrd.cc
    G4ParameterisationTubs.cc
    G4ReplicatedSlice.cc
    G4VDivisionParameterisation.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4geometrymng
    G4globman
    G4specsolids
  GLOBAL_DEPENDENCIES
    G4global
)

# List any source specific properties here
