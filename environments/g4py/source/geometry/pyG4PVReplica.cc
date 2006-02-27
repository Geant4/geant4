// $Id: pyG4PVReplica.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4PVReplica.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4PVReplica.hh"
#include "G4LogicalVolume.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4PVReplica()
{
  class_<G4PVReplica, G4PVReplica*, bases<G4VPhysicalVolume>,
    boost::noncopyable >
    ("G4PVReplica", "physical volume placement with replication", no_init)
    // constructors
    .def(init<const G4String&, G4LogicalVolume*, G4LogicalVolume*,	 
	 const EAxis, const G4int, const G4double>())
    .def(init<const G4String&, G4LogicalVolume*, G4LogicalVolume*,	 
	 const EAxis, const G4int, const G4double, const G4double>())
    .def(init<const G4String&, G4LogicalVolume*, G4VPhysicalVolume*,
	 const EAxis, const G4int, const G4double>())
    .def(init<const G4String&, G4LogicalVolume*, G4VPhysicalVolume*,
	 const EAxis, const G4int, const G4double, const G4double>())
    // ---
    .def("GetMultiplicity",  &G4PVReplica::GetMultiplicity)
    ;
}
