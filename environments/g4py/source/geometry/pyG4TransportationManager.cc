// $Id: pyG4TransportationManager.cc,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4TransportationManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"

using namespace boost::python;

// ====================================================================
// module definition
// ====================================================================
void export_G4TransportationManager()
{
  class_<G4TransportationManager, boost::noncopyable>
    ("G4TransportationManager", "manager class for transportation", no_init)
    // ---
    .def("GetTransportationManager", 
	 &G4TransportationManager::GetTransportationManager,
	 return_value_policy<reference_existing_object>())
    .staticmethod("GetTransportationManager")
    .def("GetNavigatorForTracking", 
	 &G4TransportationManager::GetNavigatorForTracking,
	 return_internal_reference<>())
    .def("GetPropagatorInField", 
	 &G4TransportationManager::GetPropagatorInField,
	 return_internal_reference<>())
    .def("GetFieldManager", 
	 &G4TransportationManager::GetFieldManager,
	 return_internal_reference<>())
    .def("SetNavigatorForTracking", 
	 &G4TransportationManager::SetNavigatorForTracking)
    .def("SetPropagatorInField", 
	 &G4TransportationManager::SetPropagatorInField)
    .def("SetFieldManager", 
	 &G4TransportationManager::SetFieldManager)
    ;
}

