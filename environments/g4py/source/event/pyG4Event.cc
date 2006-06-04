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
// $Id: pyG4Event.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
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
