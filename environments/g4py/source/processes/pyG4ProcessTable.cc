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
// $Id: pyG4ProcessTable.cc 76884 2013-11-18 12:54:03Z gcosmo $
// ====================================================================
//   pyG4ProcessTable.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
class G4UImessenger;
#include "G4ProcessTable.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4ProcessTable {

// FindProcess
G4VProcess*(G4ProcessTable::*f1_FindProcess)
  (const G4String&, const G4String&) const = &G4ProcessTable::FindProcess;

G4VProcess*(G4ProcessTable::*f2_FindProcess)
  (const G4String&, const G4ParticleDefinition*) const
  = &G4ProcessTable::FindProcess;

G4VProcess*(G4ProcessTable::*f3_FindProcess)
  (const G4String&, const G4ProcessManager*) const
  = &G4ProcessTable::FindProcess;

// FindProcesses
// raw vector pointer -> Python list conversion
list f1_FindProcesses(G4ProcessTable* procTable)
{
  list procList;
  G4ProcessVector* procVec= procTable-> FindProcesses();
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

list f2_FindProcesses(G4ProcessTable* procTable,
		      const G4ProcessManager* procManager)
{
  list procList;
  G4ProcessVector* procVec= procTable-> FindProcesses(procManager);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

list f3_FindProcesses(G4ProcessTable* procTable,
		      const G4String& pname)
{
  list procList;
  G4ProcessVector* procVec= procTable-> FindProcesses(pname);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

list f4_FindProcesses(G4ProcessTable* procTable,
		      G4ProcessType ptype)
{
  list procList;
  G4ProcessVector* procVec= procTable-> FindProcesses(ptype);
  G4int nproc= procVec-> size();
  for(G4int i=0; i< nproc; i++) {
    procList.append(&(*procVec)[i]);
  }
  return procList;
}

// SetProcessActivation
void(G4ProcessTable::*f1_SetProcessActivation)
  (const G4String&, G4bool)= &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f2_SetProcessActivation)
  (const G4String&, const G4String&, G4bool)
  = &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f3_SetProcessActivation)
  (const G4String&, G4ParticleDefinition*, G4bool)
  = &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f4_SetProcessActivation)
  (const G4String&, G4ProcessManager*, G4bool)
  = &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f5_SetProcessActivation)
  (G4ProcessType, G4bool)= &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f6_SetProcessActivation)
  (G4ProcessType, const G4String&, G4bool)
  = &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f7_SetProcessActivation)
  (G4ProcessType, G4ParticleDefinition*, G4bool)
  = &G4ProcessTable::SetProcessActivation;

void(G4ProcessTable::*f8_SetProcessActivation)
  (G4ProcessType, G4ProcessManager*, G4bool)
  = &G4ProcessTable::SetProcessActivation;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_DumpInfo, DumpInfo, 1, 2)

}

using namespace pyG4ProcessTable;

// ====================================================================
// module definition
// ====================================================================
void export_G4ProcessTable()
{
  class_<G4ProcessTable, G4ProcessTable*, boost::noncopyable>
    ("G4ProcessTable", "process table")
    // ---
    .def("GetProcessTable",  &G4ProcessTable::GetProcessTable,
         return_value_policy<reference_existing_object>())
    .staticmethod("GetProcessTable")
    .def("Length",               &G4ProcessTable::Length)
    //.def("Insert",             &G4ProcessTable::Insert)  // protected
    //.def("Remove",             &G4ProcessTable::Remove)  // protected
    // ---
    .def("FindProcess",          f1_FindProcess,
         return_value_policy<reference_existing_object>())
    .def("FindProcess",          f2_FindProcess,
         return_value_policy<reference_existing_object>())
    .def("FindProcess",          f3_FindProcess,
         return_value_policy<reference_existing_object>())
    .def("FindProcess",          f3_FindProcess,
         return_value_policy<reference_existing_object>())
    // ---
    .def("FindProcesses",        f1_FindProcesses)
    .def("FindProcesses",        f2_FindProcesses)
    .def("FindProcesses",        f3_FindProcesses)
    .def("FindProcesses",        f4_FindProcesses)
    // ---
    .def("SetProcessActivation", f1_SetProcessActivation)
    .def("SetProcessActivation", f2_SetProcessActivation)
    .def("SetProcessActivation", f3_SetProcessActivation)
    .def("SetProcessActivation", f4_SetProcessActivation)
    .def("SetProcessActivation", f5_SetProcessActivation)
    .def("SetProcessActivation", f6_SetProcessActivation)
    .def("SetProcessActivation", f7_SetProcessActivation)
    .def("SetProcessActivation", f8_SetProcessActivation)
    // ---
    .def("GetNameList",          &G4ProcessTable::GetNameList,
         return_internal_reference<>())
    .def("DumpInfo",             &G4ProcessTable::DumpInfo, f_DumpInfo())
    .def("SetVerboseLevel",      &G4ProcessTable::SetVerboseLevel)
    .def("GetVerboseLevel",      &G4ProcessTable::GetVerboseLevel)
    ;
}
