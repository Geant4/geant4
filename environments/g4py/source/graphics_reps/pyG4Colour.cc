// $Id: pyG4Colour.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Colour.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Colour.hh"
#include "G4Color.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4Colour()
{
  class_<G4Colour> ("G4Color", "color class", no_init)
    // constructors
    .def(init<>())
    .def(init<G4double>())
    .def(init<G4double, G4double>())
    .def(init<G4double, G4double, G4double>())
    .def(init<G4double, G4double, G4double, G4double>())
    .def(init<G4ThreeVector>())
    // ---
    .def("GetRed",     &G4Colour::GetRed)
    .def("GetGreen",   &G4Colour::GetGreen)
    .def("GetBlue",    &G4Colour::GetBlue)
    .def("GetAlpha",   &G4Colour::GetAlpha)
    // operators
    .def(self_ns::str(self))
    .def(self != self)
    ;

  //class_<G4Color> ("G4Color", "color class", no_init)
  //  // constructors
  //  .def(init<>())
  //  .def(init<G4double>())
  //  .def(init<G4double, G4double>())
  //  .def(init<G4double, G4double, G4double>())
  //  .def(init<G4double, G4double, G4double, G4double>())
  //  .def(init<G4ThreeVector>())
  //  // ---
  //  .def("GetRed",     &G4Colour::GetRed)
  //  .def("GetGreen",   &G4Colour::GetGreen)
  //  .def("GetBlue",    &G4Colour::GetBlue)
  //  .def("GetAlpha",   &G4Colour::GetAlpha)
  //  // operators
  //  .def(self_ns::str(self))
  //  .def(self != self)
  //  ;
}

