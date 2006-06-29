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
// $Id: pyG4Trd.cc,v 1.4 2006-06-29 15:32:42 gunter Exp $
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

