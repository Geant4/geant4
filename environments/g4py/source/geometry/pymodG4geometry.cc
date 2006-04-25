// $Id: pymodG4geometry.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pymodG4geometry.cc [Geant4Py module]
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4VTouchable();
void export_G4TouchableHistory();
void export_G4VPhysicalVolume();
void export_G4PVPlacement();
void export_G4PVReplica();
void export_G4LogicalVolume();
void export_G4Region();
void export_G4VSolid();
void export_G4Box();
void export_G4Cons();
void export_G4Para();
void export_G4Torus();
void export_G4Trd();
void export_G4Orb();
void export_G4Sphere();
void export_G4Trap();
void export_G4Tubs();
void export_G4TransportationManager();
void export_G4FieldManager();
void export_G4Field();
void export_G4MagneticField();
void export_G4UniformMagField();
void export_G4ChordFinder();

BOOST_PYTHON_MODULE(G4geometry)
{
  export_G4VTouchable();
  export_G4TouchableHistory();
  export_G4VPhysicalVolume();
  export_G4PVPlacement();
  export_G4PVReplica();
  export_G4LogicalVolume();
  export_G4Region();
  export_G4VSolid();
  export_G4Box();
  export_G4Cons();
  export_G4Para();
  export_G4Torus();
  export_G4Trd();
  export_G4Orb();
  export_G4Sphere();
  export_G4Trap();
  export_G4Tubs();
  export_G4TransportationManager();
  export_G4FieldManager();
  export_G4Field();
  export_G4MagneticField();
  export_G4UniformMagField();
  export_G4ChordFinder();
}

