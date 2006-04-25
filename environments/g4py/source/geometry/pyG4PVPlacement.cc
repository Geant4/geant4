// $Id: pyG4PVPlacement.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4PVPlacement.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4PVPlacement()
{
  class_<G4PVPlacement, G4PVPlacement*, bases<G4VPhysicalVolume>, 
    boost::noncopyable >
    ("G4PVPlacement", "physical volume placement", no_init)
    // ---
    .def(init<G4RotationMatrix*, const G4ThreeVector&,
	 G4LogicalVolume*, const G4String&,
	 G4LogicalVolume*, G4bool, G4int>())
    .def(init<const G4Transform3D&, G4LogicalVolume*,
	 const G4String&, G4LogicalVolume*, G4bool, G4int>())
    .def(init<G4RotationMatrix*, const G4ThreeVector&,
	 const G4String, G4LogicalVolume*,
	 G4VPhysicalVolume*, G4bool, G4int>())
    .def(init<const G4Transform3D&, const G4String&,
	 G4LogicalVolume*, G4VPhysicalVolume*, G4bool, G4int>())
    ;
}

