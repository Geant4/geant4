// $Id: pyG4Box.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Box.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Box.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Box()
{
  class_<G4Box, G4Box*, bases<G4VSolid> >
    ("G4Box", "box solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double>())
    // ---
    .def("GetXHalfLength",   &G4Box::GetXHalfLength)
    .def("GetYHalfLength",   &G4Box::GetYHalfLength)
    .def("GetZHalfLength",   &G4Box::GetZHalfLength)
    .def("SetXHalfLength",   &G4Box::SetXHalfLength)
    .def("SetYHalfLength",   &G4Box::SetYHalfLength)
    .def("SetZHalfLength",   &G4Box::SetZHalfLength)
    .def("GetCubicVolume",   &G4Box::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

