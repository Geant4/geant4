# - G4geomdivision module build definition

# Define the Geant4 Module.
geant4_add_module(G4geomdivision
  PUBLIC_HEADERS
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
    G4VDivisionParameterisation.cc)

geant4_module_link_libraries(G4geomdivision
  PUBLIC G4specsolids G4volumes G4globman G4geometrymng G4hepgeometry
  PRIVATE G4csg)
