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
// $Id: pyG4VProcess.cc,v 1.4 2006-06-29 15:35:03 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VProcess.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VProcess.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VProcess {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_StorePhysicsTable, 
				       StorePhysicsTable, 2, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_RetrievePhysicsTable, 
				       RetrievePhysicsTable, 2, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetPhysicsTableFileName, 
				       GetPhysicsTableFileName, 3, 4);
};

using namespace pyG4VProcess;

// ====================================================================
// module definition
// ====================================================================
void export_G4VProcess()
{
  class_<G4VProcess, G4VProcess*, boost::noncopyable>
    ("G4VProcess", "base class for process", no_init)
    // ---
    // Note that only limited methods are exposed.
    .def("SetPILfactor",         &G4VProcess::SetPILfactor)
    .def("GetPILfactor",         &G4VProcess::GetPILfactor)
    .def("IsApplicable",         &G4VProcess::IsApplicable)
    .def("BuildPhysicsTable",    &G4VProcess::BuildPhysicsTable)
    .def("PreparePhysicsTable",  &G4VProcess::PreparePhysicsTable)
    .def("StorePhysicsTable",    &G4VProcess::StorePhysicsTable,
	 f_StorePhysicsTable())
    .def("RetrievePhysicsTable", &G4VProcess::RetrievePhysicsTable,
	 f_RetrievePhysicsTable())
    .def("GetPhysicsTableFileName", &G4VProcess::GetPhysicsTableFileName,
	 f_GetPhysicsTableFileName()
	 [return_value_policy<return_by_value>()])
    .def("GetProcessName",       &G4VProcess::GetProcessName,
	 return_value_policy<return_by_value>())  
    .def("GetProcessType",       &G4VProcess::GetProcessType)    
    .def("DumpInfo",             &G4VProcess::DumpInfo)
    .def("SetVerboseLevel",      &G4VProcess::SetVerboseLevel)
    .def("GetVerboseLevel",      &G4VProcess::GetVerboseLevel)
    ;
}
