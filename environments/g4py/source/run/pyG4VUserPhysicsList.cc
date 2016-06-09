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
// $Id: pyG4VUserPhysicsList.cc,v 1.5 2006-06-29 15:35:30 gunter Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4VUserPhysicsList.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VUserPhysicsList.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4VUserPhysicsList {

struct CB_G4VUserPhysicsList :
  G4VUserPhysicsList, wrapper<G4VUserPhysicsList> {

  void ConstructParticle() {
    get_override("ConstructParticle")();
  }

  void ConstructProcess() {
    get_override("ConstructProcess")();
  }

  void SetCuts() {
    get_override("SetCuts")();
  }
};

// SetCutValue
void (G4VUserPhysicsList::*f1_SetCutValue)(G4double, const G4String&)
  = &G4VUserPhysicsList::SetCutValue;
void (G4VUserPhysicsList::*f2_SetCutValue)(G4double, const G4String&,
					   const G4String&)
  = &G4VUserPhysicsList::SetCutValue;

// SetParticleCuts
void (G4VUserPhysicsList::*f1_SetParticleCuts)(G4double,
                                               G4ParticleDefinition*,
                                               G4Region*)
  = &G4VUserPhysicsList::SetParticleCuts;
void (G4VUserPhysicsList::*f2_SetParticleCuts)(G4double,
                                               G4ParticleDefinition*,
                                               G4Region*)
  = &G4VUserPhysicsList::SetParticleCuts;

// StorePhysicsTable
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_StorePhysicsTable,
				       StorePhysicsTable, 0, 1);
// SetParticleCuts
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_SetParticleCuts,
				       SetParticleCuts, 2, 3);

}

using namespace pyG4VUserPhysicsList;

// ====================================================================
// module definition
// ====================================================================
void export_G4VUserPhysicsList()
{
  class_<CB_G4VUserPhysicsList, CB_G4VUserPhysicsList*, boost::noncopyable>
    ("G4VUserPhysicsList", "base class of user physics list")
    // ---
    .def("ConstructParticle",
	 pure_virtual(&G4VUserPhysicsList::ConstructParticle))
    .def("ConstructProcess",
	 pure_virtual(&G4VUserPhysicsList::ConstructProcess))
    .def("SetCuts",
	 pure_virtual(&G4VUserPhysicsList::SetCuts))
    // ---
    .def("SetDefaultCutValue",    &G4VUserPhysicsList::SetDefaultCutValue)
    .def("GetDefaultCutValue",    &G4VUserPhysicsList::GetDefaultCutValue)
    // ---
    .def("StorePhysicsTable",     &G4VUserPhysicsList::StorePhysicsTable,
	 f_StorePhysicsTable())
    .def("IsPhysicsTableRetrieved", 
         &G4VUserPhysicsList::IsPhysicsTableRetrieved)
    .def("IsStoredInAscii",       &G4VUserPhysicsList::IsStoredInAscii)
    .def("GetPhysicsTableDirectory", 
         &G4VUserPhysicsList::GetPhysicsTableDirectory,
         return_value_policy<return_by_value>())
    .def("SetStoredInAscii",      &G4VUserPhysicsList::SetStoredInAscii)
    .def("ResetStoredInAscii",    &G4VUserPhysicsList::ResetStoredInAscii)
    // ---
    .def("DumpList",              &G4VUserPhysicsList::DumpList)

    .def("DumpCutValuesTable",    &G4VUserPhysicsList::DumpCutValuesTable)
    .def("DumpCutValuesTableIfRequested", 
         &G4VUserPhysicsList::DumpCutValuesTableIfRequested)
    .def("SetCutValue",           f1_SetCutValue)
    .def("SetCutValue",           f2_SetCutValue)
    .def("SetParticleCuts",       f1_SetParticleCuts, f_SetParticleCuts())
    .def("SetParticleCuts",       f2_SetParticleCuts, f_SetParticleCuts())
    // ---
    .def("SetVerboseLevel",       &G4VUserPhysicsList::SetVerboseLevel)
    .def("GetVerboseLevel",       &G4VUserPhysicsList::GetVerboseLevel)
    .def("SetCutsWithDefault",    &G4VUserPhysicsList::SetCutsWithDefault)
    .def("SetCutsForRegion",      &G4VUserPhysicsList::SetCutsForRegion)
    .def("GetApplyCuts",          &G4VUserPhysicsList::GetApplyCuts)
    ;

  // Note that exposed items are limited,
  // because this class object is mainly for internal uses or obsolete.
  // Construct
  // BuildPhysicsTable
  // PreparePhysicsTable
  // SetPhysicsTableRetrieved
  // ReSetPhysicsTableRetrieved
  // SetApplyCuts
  // DumpCutValues (obsolete)
  // ResetCuts;
}
