# - G4hepgeometry module build definition

geant4_add_module(G4hepgeometry
  PUBLIC_HEADERS
    G4LorentzRotation.hh
    G4LorentzVector.hh
    G4Normal3D.hh
    G4Plane3D.hh
    G4Point3D.hh
    G4Transform3D.hh
    G4Vector3D.hh
    geomdefs.hh)

geant4_module_link_libraries(G4hepgeometry INTERFACE G4globman)
