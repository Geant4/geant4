// $Id: pyG4Trd.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Trd.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Trd.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Trd()
{
  class_<G4Trd, G4Trd*, bases<G4VSolid> >
    ("G4Trd", "Trapezoild solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double>())
    // ---
    .def("GetXHalfLength1", &G4Trd::GetXHalfLength1)
    .def("GetXHalfLength2", &G4Trd::GetXHalfLength2)
    .def("GetYHalfLength1", &G4Trd::GetYHalfLength1)
    .def("GetYHalfLength2", &G4Trd::GetYHalfLength2)
    .def("GetZHalfLength",  &G4Trd::GetZHalfLength)
    .def("SetXHalfLength1", &G4Trd::SetXHalfLength1)
    .def("SetXHalfLength2", &G4Trd::SetXHalfLength2)
    .def("SetYHalfLength1", &G4Trd::SetYHalfLength1)
    .def("SetYHalfLength2", &G4Trd::SetYHalfLength2)
    .def("SetZHalfLength",  &G4Trd::SetZHalfLength)
    .def("GetCubicVolume",  &G4Trd::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

