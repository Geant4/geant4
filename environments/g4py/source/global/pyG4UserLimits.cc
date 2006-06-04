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
// $Id: pyG4UserLimits.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
// $Name: not supported by cvs2svn $
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
