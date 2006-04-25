// $Id: pyG4Trap.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Trap.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Trap.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Trap()
{
  class_<G4Trap, G4Trap*, bases<G4VSolid> >
    ("G4Trap", "Generic trapezoild soild class", no_init)
    // constructors
    .def(init<const G4String&>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double, G4double, G4double, G4double,
	 G4double, G4double, G4double>())
    // ---
    .def("GetZHalfLength",   &G4Trap::GetZHalfLength)
    .def("GetYHalfLength1",  &G4Trap::GetYHalfLength1)
    .def("GetXHalfLength1",  &G4Trap::GetXHalfLength1)
    .def("GetXHalfLength2",  &G4Trap::GetXHalfLength2)
    .def("GetTanAlpha1",     &G4Trap::GetTanAlpha1)
    .def("GetYHalfLength2",  &G4Trap::GetYHalfLength2)
    .def("GetXHalfLength3",  &G4Trap::GetXHalfLength3)
    .def("GetXHalfLength4",  &G4Trap::GetXHalfLength4)
    .def("GetTanAlpha2",     &G4Trap::GetTanAlpha2)
    .def("GetSidePlane",     &G4Trap::GetSidePlane)
    .def("GetSymAxis",       &G4Trap::GetSymAxis)
    .def("SetAllParameters", &G4Trap::SetAllParameters)
    .def("GetCubicVolume",   &G4Trap::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
