//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: pyG4Colour.cc,v 1.4 2006-06-29 15:33:54 gunter Exp $
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

