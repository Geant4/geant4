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
// $Id: pyG4ProcessManager.cc,v 1.4 2006-06-29 15:34:55 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ProcessManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4ProcessManager.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ProcessManager {

// GetProcessList()
// raw vector pointer -> Python list conversion
list f_GetProcessList(const G4ProcessManager* procMgr)
{
  list procList;
  G4ProcessVector* procVec= procMgr-> GetProcessList();
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

// GetProcessVector()
list f_GetProcessVector(const G4ProcessManager* procMgr, 
			G4ProcessVectorDoItIndex idx,
			G4ProcessVectorTypeIndex typ= typeGPIL )
{
  list procList;
  G4ProcessVector* procVec= procMgr-> GetProcessVector(idx, typ);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(g_GetProcessVector,
				f_GetProcessVector, 2, 3);

// GetAtRestProcessVector()
list f_GetAtRestProcessVector(const G4ProcessManager* procMgr, 
			      G4ProcessVectorTypeIndex typ= typeGPIL )
{
  list procList;
  G4ProcessVector* procVec= procMgr-> GetAtRestProcessVector(typ);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(g_GetAtRestProcessVector,
				f_GetAtRestProcessVector, 1, 2);

// GetAlongStepProcessVector()
list f_GetAlongStepProcessVector(const G4ProcessManager* procMgr, 
				 G4ProcessVectorTypeIndex typ= typeGPIL )
{
  list procList;
  G4ProcessVector* procVec= procMgr-> GetAlongStepProcessVector(typ);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(g_GetAlongStepProcessVector,
				f_GetAlongStepProcessVector, 1, 2);

// GetPostStepProcessVector()
list f_GetPostStepProcessVector(const G4ProcessManager* procMgr, 
				G4ProcessVectorTypeIndex typ= typeGPIL )
{
  list procList;
  G4ProcessVector* procVec= procMgr-> GetPostStepProcessVector(typ);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(g_GetPostStepProcessVector,
				f_GetPostStepProcessVector, 1, 2);


// GetProcessVectorIndex...
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetProcessVectorIndex,
				       GetProcessVectorIndex, 2, 3);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetAtRestIndex,
				       GetAtRestIndex, 1, 2);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetAlongStepIndex,
				       GetAlongStepIndex, 1, 2);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_GetPostStepIndex,
				       GetPostStepIndex, 1, 2);
// AddProcess...
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_AddProcess, AddProcess, 1, 4);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_AddRestProcess, 
				       AddRestProcess, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_AddDiscreteProcess, 
				       AddDiscreteProcess, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_AddContinuousProcess, 
				       AddContinuousProcess, 1, 2);
// SetProcessOrdering
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetProcessOrdering, 
				       SetProcessOrdering, 2, 3);

// RemoveProcess...
G4VProcess*(G4ProcessManager::*f1_RemoveProcess)(G4VProcess*)
  = &G4ProcessManager::RemoveProcess;

G4VProcess*(G4ProcessManager::*f2_RemoveProcess)(G4int)
  = &G4ProcessManager::RemoveProcess;

// Set/GetProcessActivation...
G4VProcess*(G4ProcessManager::*f1_SetProcessActivation)(G4VProcess*, G4bool)
  = &G4ProcessManager::SetProcessActivation;

G4VProcess*(G4ProcessManager::*f2_SetProcessActivation)(G4int, G4bool)
  = &G4ProcessManager::SetProcessActivation;

G4bool(G4ProcessManager::*f1_GetProcessActivation)(G4VProcess*) const
  = &G4ProcessManager::GetProcessActivation;

G4bool(G4ProcessManager::*f2_GetProcessActivation)(G4int) const
  = &G4ProcessManager::GetProcessActivation;

};

using namespace pyG4ProcessManager;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProcessManager()
{
  class_<G4ProcessManager, G4ProcessManager*, boost::noncopyable>
    ("G4ProcessManager", "process manager class", no_init)    
    // ---
    .def("GetProcessList",            f_GetProcessList)
    .def("GetProcessListLength",      &G4ProcessManager::GetProcessListLength)
    .def("GetProcessIndex",           &G4ProcessManager::GetProcessIndex)
    .def("GetProcessVector",          f_GetProcessVector, 
	                              g_GetProcessVector())
    .def("GetAtRestProcessVector",    f_GetAtRestProcessVector,
	                              g_GetAtRestProcessVector())
    .def("GetAlongStepProcessVector", f_GetAlongStepProcessVector,
	                              g_GetAlongStepProcessVector())
    .def("GetPostStepProcessVector",  f_GetPostStepProcessVector,
	                              g_GetPostStepProcessVector())
    .def("GetProcessVectorIndex",     
	 &G4ProcessManager::GetProcessVectorIndex,
	 f_GetProcessVectorIndex())
    .def("GetAtRestIndex",            &G4ProcessManager::GetAtRestIndex,
	                              f_GetAtRestIndex())
    .def("GetAlongStepIndex",         &G4ProcessManager::GetAlongStepIndex,
	                              f_GetAlongStepIndex())
    .def("GetPostStepIndex",          &G4ProcessManager::GetPostStepIndex,
	                              f_GetPostStepIndex())
    // ----
    .def("AddProcess",            &G4ProcessManager::AddProcess,
	                          f_AddProcess())
    .def("AddRestProcess",        &G4ProcessManager::AddRestProcess,
	                          f_AddRestProcess())
    .def("AddDiscreteProcess",    &G4ProcessManager::AddDiscreteProcess,
	                          f_AddDiscreteProcess())
    .def("AddContinuousProcess",  &G4ProcessManager::AddContinuousProcess,
	                          f_AddContinuousProcess())
    // ---
    .def("GetProcessOrdering",    &G4ProcessManager::GetProcessOrdering)
    .def("SetProcessOrdering",    &G4ProcessManager::SetProcessOrdering,
             	                  f_SetProcessOrdering()) 
    .def("SetProcessOrderingToFirst", 
	 &G4ProcessManager::SetProcessOrderingToFirst)
    .def("SetProcessOrderingToLast", 
	 &G4ProcessManager::SetProcessOrderingToLast)
    // ---
    .def("RemoveProcess",         f1_RemoveProcess,
         return_value_policy<reference_existing_object>())
    .def("RemoveProcess",         f2_RemoveProcess,
         return_value_policy<reference_existing_object>())
    // ---
    .def("SetProcessActivation",  f1_SetProcessActivation,
	 return_value_policy<reference_existing_object>())
    .def("SetProcessActivation",  f2_SetProcessActivation,
	 return_value_policy<reference_existing_object>())
    .def("GetProcessActivation",  f1_GetProcessActivation)
    .def("GetProcessActivation",  f2_GetProcessActivation)
    // ---
    .def("GetParticleType",       &G4ProcessManager::GetParticleType,
         return_internal_reference<>())
    .def("SetParticleType",       &G4ProcessManager::SetParticleType)
    .def("DumpInfo",              &G4ProcessManager::DumpInfo)
    .def("SetVerboseLevel",       &G4ProcessManager::SetVerboseLevel)
    .def("GetVerboseLevel",       &G4ProcessManager::GetVerboseLevel)
    ;

  // enums...
  enum_<G4ProcessVectorTypeIndex>("G4ProcessVectorTypeIndex")
    .value("typeGPIL", typeGPIL)
    .value("typeGPIL", typeDoIt)
    ;

  enum_<G4ProcessVectorDoItIndex>("G4ProcessVectorDoItIndex")
    .value("idxAll",       idxAll)
    .value("idxAtRest",    idxAtRest)
    .value("idxAlongStep", idxAlongStep)
    .value("idxPostStep",  idxPostStep)
    ;

  enum_<G4ProcessVectorOrdering>("G4ProcessVectorOrdering")
    .value("ordInActive", ordInActive)
    .value("ordDefault",  ordDefault)
    .value("ordLast",     ordLast)
    ;
}
