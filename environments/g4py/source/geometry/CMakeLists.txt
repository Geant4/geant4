# - build library

# library
geant4_add_pymodule(pyG4geometry
  pyG4BooleanSolid.cc
  pyG4Box.cc
  pyG4ChordFinder.cc
  pyG4Cons.cc
  pyG4Ellipsoid.cc
  pyG4EllipticalCone.cc
  pyG4EllipticalTube.cc
  pyG4Field.cc
  pyG4FieldManager.cc
  pyG4GeometryManager.cc
  pyG4Hype.cc
  pyG4IntersctionSolid.cc
  pyG4LogicalVolume.cc
  pyG4MagneticField.cc
  pyG4Navigator.cc
  pyG4Orb.cc
  pyG4PVPlacement.cc
  pyG4PVReplica.cc
  pyG4Para.cc
  pyG4Polycone.cc
  pyG4Polyhedra.cc
  pyG4Region.cc
  pyG4Sphere.cc
  pyG4SubtractionSolid.cc
  pyG4Tet.cc
  pyG4Torus.cc
  pyG4TouchableHistory.cc
  pyG4TransportationManager.cc
  pyG4Trap.cc
  pyG4Trd.cc
  pyG4Tubs.cc
  pyG4TwistedBox.cc
  pyG4TwistedTrap.cc
  pyG4TwistedTrd.cc
  pyG4TwistedTubs.cc
  pyG4UniformMagField.cc
  pyG4UnionSolid.cc
  pyG4VPhysicalVolume.cc
  pyG4VSolid.cc
  pyG4VTouchable.cc
  pymodG4geometry.cc
)

target_link_libraries(pyG4geometry PRIVATE G4geometry G4track G4processes)
install(TARGETS pyG4geometry DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
