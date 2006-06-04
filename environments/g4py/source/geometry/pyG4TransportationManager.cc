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
// $Id: pyG4TransportationManager.cc,v 1.3 2006-06-04 21:34:29 kmura Exp $
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

