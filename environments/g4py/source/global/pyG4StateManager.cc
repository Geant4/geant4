// $Id: pyG4StateManager.cc,v 1.2 2006-04-25 08:09:45 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4StateManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4StateManager.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4StateManager()
{
  class_<G4StateManager, boost::noncopyable>
    ("G4StateManager", "state manager", no_init)
    .def("GetStateManager",    &G4StateManager::GetStateManager,
	 "Get an instance of G4StateManager",
         return_value_policy<reference_existing_object>())
    .staticmethod("GetStateManager")
    // ---
    .def("GetCurrentState",    &G4StateManager::GetCurrentState)
    .def("GetPreviousState",   &G4StateManager::GetPreviousState)
    .def("GetStateString",     &G4StateManager::GetStateString)
    ;
}

