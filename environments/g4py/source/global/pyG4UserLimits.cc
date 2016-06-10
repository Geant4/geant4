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
// $Id: pyG4UserLimits.cc 66892 2013-01-17 10:57:59Z gunter $
// ====================================================================
//   pyG4UserLimits.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4UserLimits.hh"
#include "G4Track.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4UserLimits()
{
  class_<G4UserLimits, G4UserLimits*>
    ("G4UserLimits", "user step limitations")
    // ---
    .def(init<G4double>())
    .def(init<G4double, G4double>())
    .def(init<G4double, G4double, G4double>())
    .def(init<G4double, G4double, G4double, G4double>())
    .def(init<G4double, G4double, G4double, G4double, G4double>())
    // ---
    .def(init<const G4String&>())
    .def(init<const G4String&, G4double>())
    .def(init<const G4String&, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, G4double>())
    .def(init<const G4String&, G4double, G4double, G4double, 
	 G4double, G4double>())
    // ---
    .def("GetUserMaxTrackLength", &G4UserLimits::GetUserMaxTrackLength)
    .def("GetUserMaxTime",        &G4UserLimits::GetUserMaxTime)
    .def("GetUserMinEkine",       &G4UserLimits::GetUserMinEkine)
    .def("GetUserMinRange",       &G4UserLimits::GetUserMinRange)
    // ---
    .def("SetMaxAllowedStep",     &G4UserLimits::SetMaxAllowedStep)
    .def("SetUserMaxTrackLength", &G4UserLimits::SetUserMaxTrackLength)
    .def("SetUserMaxTime",        &G4UserLimits::SetUserMaxTime)
    .def("SetUserMinEkine",       &G4UserLimits::SetUserMinEkine)
    .def("SetUserMinRange",       &G4UserLimits::SetUserMinRange)
    // ---
    .def("GetType",               &G4UserLimits::GetType,
         return_internal_reference<>())
    .def("SetType",               &G4UserLimits::SetType)
    ;
}
