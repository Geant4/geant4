#------------------------------------------------------------------------------
# Module : G4hepgeometry
# Package: Geant4.src.G4global.G4hepgeometry
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4hepgeometry
  HEADERS
    G4LorentzRotation.hh
    G4LorentzVector.hh
    G4Normal3D.hh
    G4Plane3D.hh
    G4Point3D.hh
    G4Transform3D.hh
    G4Vector3D.hh
    geomdefs.hh
  GRANULAR_DEPENDENCIES
    G4globman
)

# List any source specific properties here
