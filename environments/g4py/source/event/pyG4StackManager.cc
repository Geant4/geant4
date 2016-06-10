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

// $Id: pyG4StackManager.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4StackManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4StackManager.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4StackManager {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_ClearWaitingStack,
                                       ClearWaitingStack, 0, 1)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetNWaitingTrack,
                                       GetNWaitingTrack, 0, 1)
}

using namespace pyG4StackManager;

// ====================================================================
// module definition
// ====================================================================
void export_G4StackManager()
{
  class_<G4StackManager, boost::noncopyable>
    ("G4StackManager", "stack management class")
    // ---
    // Note that exposed items are limited, because this class object
    // is mainly for internal uses.
    .def("ReClassify",          &G4StackManager::ReClassify)
    .def("clear",               &G4StackManager::clear)
    .def("ClearUrgentStack",    &G4StackManager::ClearUrgentStack)
    .def("ClearWaitingStack",   &G4StackManager::ClearWaitingStack,
	 f_ClearWaitingStack())
    .def("ClearPostponeStack",  &G4StackManager::ClearPostponeStack)
    .def("GetNTotalTrack",      &G4StackManager::GetNTotalTrack)
    .def("GetNUrgentTrack",     &G4StackManager::GetNUrgentTrack)
    .def("GetNWaitingTrack",    &G4StackManager::GetNWaitingTrack,
	 f_GetNWaitingTrack())
    .def("SetVerboseLevel",     &G4StackManager::SetVerboseLevel)
    ;
}
