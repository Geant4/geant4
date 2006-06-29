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
// $Id: pyG4Event.cc,v 1.4 2006-06-29 15:31:29 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Event.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Event.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Event {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetPrimaryVertex, 
				       GetPrimaryVertex, 0, 1);
};

using namespace pyG4Event;

// ====================================================================
// module definition
// ====================================================================
void export_G4Event()
{
  class_<G4Event, G4Event*>("G4Event", "event class")
    .def(init<G4int>())
    // ---
    .def("Print",              &G4Event::Print)
    .def("Draw",               &G4Event::Draw)
    .def("SetEventID",         &G4Event::SetEventID)
    .def("GetEventID",         &G4Event::GetEventID)
    .def("SetEventAborted",    &G4Event::SetEventAborted)
    .def("IsAborted",          &G4Event::IsAborted)
    // ---
    .def("AddPrimaryVertex",   &G4Event::AddPrimaryVertex)
    .def("GetNumberOfPrimaryVertex", &G4Event::GetNumberOfPrimaryVertex)
    .def("GetPrimaryVertex",   &G4Event::GetPrimaryVertex,
	 f_GetPrimaryVertex()[return_internal_reference<>()])
    // ---
    .def("GetTrajectoryContainer", &G4Event::GetTrajectoryContainer,
	 return_internal_reference<>())
    .def("SetUserInformation", &G4Event::SetUserInformation)
    .def("GetUserInformation", &G4Event::GetUserInformation,
	 return_internal_reference<>())
    ;

  // reduced functionality...
  //.def("SetHCofThisEvent",   &G4Event::SetHCofThisEvent)
  //.def("GetHCofThisEvent",   &G4Event::SetHCofThisEvent,
  //     return_internal_reference<>())
  //.def("SetDCofThisEvent",   &G4Event::SetHCofThisEvent)
  //.def("GetDCofThisEvent",   &G4Event::SetHCofThisEvent,
  //     return_internal_reference<>())

}
