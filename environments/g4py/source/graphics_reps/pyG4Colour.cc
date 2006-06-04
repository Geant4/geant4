//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: pyG4Colour.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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

