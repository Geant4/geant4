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
