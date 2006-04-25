// $Id: pyG4Orb.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Orb.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Orb.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Orb()
{
  class_<G4Orb, G4Orb*, bases<G4VSolid> >
    ("G4Orb", "Orb solid class", no_init)
    // constructors
    .def(init<const G4String&, G4double>())
    // ---
    .def("GetRadius", &G4Orb::GetRadius)
    .def("SetRadius", &G4Orb::SetRadius)
    .def("GetCubicVolume", &G4Orb::GetCubicVolume)
    // operators
    .def(self_ns::str(self))
    ;
}
