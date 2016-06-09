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
// $Id: pyG4Trap.cc,v 1.4 2006/06/29 15:32:39 gunter Exp $
// $Name: geant4-08-01 $
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
