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

// $Id: pyG4StackManager.cc,v 1.3 2006-06-04 21:34:28 kmura Exp $
// $Name: not supported by cvs2svn $
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
                                       ClearWaitingStack, 0, 1);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetNWaitingTrack,
                                       GetNWaitingTrack, 0, 1);
};

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
