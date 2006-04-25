// $Id: pyG4Para.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Para.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Para.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Para()
{
  class_<G4Para, G4Para*, bases<G4VSolid> >
    ("G4Para", "Skewed box sold class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
	 G4double, G4double, G4double>())
    // ---
    .def("GetZHalfLength",    &G4Para::GetZHalfLength)
    .def("GetSymAxis",        &G4Para::GetSymAxis)
    .def("GetYHalfLength",    &G4Para::GetYHalfLength)
    .def("GetXHalfLength",    &G4Para::GetXHalfLength)
    .def("GetTanAlpha",       &G4Para::GetTanAlpha)
    .def("SetXHalfLength",    &G4Para::SetXHalfLength)
    .def("SetYHalfLength",    &G4Para::SetYHalfLength)
    .def("SetZHalfLength",    &G4Para::SetZHalfLength)
    .def("SetAlpha",          &G4Para::SetAlpha)
    .def("SetTanAlpha",       &G4Para::SetTanAlpha)
    .def("SetThetaAndPhi",    &G4Para::SetThetaAndPhi)
    .def("SetAllParameters",  &G4Para::SetAllParameters)
    .def("GetCubicVolume",    &G4Para::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}

