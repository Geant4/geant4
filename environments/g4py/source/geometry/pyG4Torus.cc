// $Id: pyG4Torus.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Torus.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Torus.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Torus()
{
  class_<G4Torus, G4Torus*, bases<G4VSolid> >
    ("G4Torus", "Torus solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double, G4double, G4double,
         G4double, G4double>())
    // ---
    .def("GetRmin", &G4Torus::GetRmin)
    .def("GetRmax", &G4Torus::GetRmax)
    .def("GetRtor", &G4Torus::GetRtor)
    .def("GetSPhi", &G4Torus::GetSPhi)
    .def("GetDPhi", &G4Torus::GetDPhi)
    .def("GetCubicVolume", &G4Torus::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
