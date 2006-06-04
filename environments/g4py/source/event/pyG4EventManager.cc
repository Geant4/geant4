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
// $Id: pyG4EventManager.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4EventManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4EventManager.hh"
#include "G4Event.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4EventManager()
{
  class_<G4EventManager, boost::noncopyable>
    ("G4EventManager", "event manager class")
    .def("GetEventManager",  &G4EventManager::GetEventManager,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetEventManager")
    // ---
    .def("GetConstCurrentEvent", &G4EventManager::GetConstCurrentEvent,
         return_internal_reference<>())
    .def("GetNonconstCurrentEvent", 
	 &G4EventManager::GetNonconstCurrentEvent,
         return_internal_reference<>())
    .def("AbortCurrentEvent",    &G4EventManager::AbortCurrentEvent)
    .def("SetNumberOfAdditionalWaitingStacks", 
	 &G4EventManager::SetNumberOfAdditionalWaitingStacks)
    .def("GetStackManager",      &G4EventManager::GetStackManager,
	 return_value_policy<reference_existing_object>())
    .def("GetTrackingManager",   &G4EventManager::GetTrackingManager,
	 return_value_policy<reference_existing_object>())
    .def("GetVerboseLevel",      &G4EventManager::GetVerboseLevel)
    .def("SetVerboseLevel",      &G4EventManager::SetVerboseLevel)
    .def("SetUserInformation",   &G4EventManager::SetUserInformation)
    .def("GetUserInformation",   &G4EventManager::GetUserInformation,
         return_value_policy<reference_existing_object>())
    ;

  // Note that exposed items are limited, 
  // because this class object is mainly for internal uses.
  // ProcessOneEvent
  // SetUserAction
  // GetUserXXXAction
  // GetPrimaryTransformer
  // SetPrimaryTransformer

}
