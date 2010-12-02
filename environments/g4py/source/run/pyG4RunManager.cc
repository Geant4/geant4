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
// $Id: pyG4RunManager.cc,v 1.6 2010-12-02 08:23:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4RunManager.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4Version.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4Region.hh"
#include "G4Run.hh"
#include "G4Event.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4RunManager {

// SetUserInitialization()
void (G4RunManager::*f1_SetUserInitialization)(G4VUserDetectorConstruction*)
  = &G4RunManager::SetUserInitialization;
void (G4RunManager::*f2_SetUserInitialization)(G4VUserPhysicsList*)
  = &G4RunManager::SetUserInitialization;

// SetUserAction()
void (G4RunManager::*f1_SetUserAction)(G4UserRunAction*)
  = &G4RunManager::SetUserAction;
void (G4RunManager::*f2_SetUserAction)(G4VUserPrimaryGeneratorAction*)
  = &G4RunManager::SetUserAction;
void (G4RunManager::*f3_SetUserAction)(G4UserEventAction*)
  = &G4RunManager::SetUserAction;
void (G4RunManager::*f4_SetUserAction)(G4UserStackingAction*)
  = &G4RunManager::SetUserAction;
void (G4RunManager::*f5_SetUserAction)(G4UserTrackingAction*)
  = &G4RunManager::SetUserAction;
void (G4RunManager::*f6_SetUserAction)(G4UserSteppingAction*)
  = &G4RunManager::SetUserAction;

// DumpRegion
#if G4VERSION_NUMBER >= 932
void (G4RunManager::*f1_DumpRegion)(const G4String&) const
  = &G4RunManager::DumpRegion;
#else
void (G4RunManager::*f1_DumpRegion)(G4String) const
  = &G4RunManager::DumpRegion;
#endif
void (G4RunManager::*f2_DumpRegion)(G4Region*) const
  = &G4RunManager::DumpRegion;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_DumpRegion, DumpRegion, 0, 1);

// BeamOn()
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_BeamOn, BeamOn, 1, 3);

// AbortRun()
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_AbortRun, AbortRun, 0, 1);

// DefineWorldVolume()
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f_DefineWorldVolume,
               DefineWorldVolume, 1, 2);

};

using namespace pyG4RunManager;

// ====================================================================
// module definition
// ====================================================================
void export_G4RunManager()
{
  class_<G4RunManager>("G4RunManager", "run manager class")
    // ---
    .def("GetRunManager", &G4RunManager::GetRunManager,
   "Get an instance of G4RunManager",
   return_value_policy<reference_existing_object>())
    .staticmethod("GetRunManager")
    // ---
    .def("SetVerboseLevel", &G4RunManager::SetVerboseLevel)
    .def("GetVerboseLevel", &G4RunManager::GetVerboseLevel)
    // ---
    .def("Initialize",      &G4RunManager::Initialize)
    .def("BeamOn",          &G4RunManager::BeamOn,
   f_BeamOn((arg("n_event"), arg("macroFile")=0,
       arg("n_select")=-1),
      "Starts event loop."))
    // ---
    .def("SetUserInitialization", f1_SetUserInitialization)
    .def("SetUserInitialization", f2_SetUserInitialization)
    .def("SetUserAction",         f1_SetUserAction)
    .def("SetUserAction",         f2_SetUserAction)
    .def("SetUserAction",         f3_SetUserAction)
    .def("SetUserAction",         f4_SetUserAction)
    .def("SetUserAction",         f5_SetUserAction)
    .def("SetUserAction",         f6_SetUserAction)
    // ---
    .def("GetUserDetectorConstruction",
   &G4RunManager::GetUserDetectorConstruction,
   return_internal_reference<>())
    .def("GetUserPhysicsList",
   &G4RunManager::GetUserPhysicsList,
   return_internal_reference<>())
    .def("GetUserPrimaryGeneratorAction",
   &G4RunManager::GetUserPrimaryGeneratorAction,
   return_internal_reference<>())
    .def("GetUserRunAction",      &G4RunManager::GetUserRunAction,
   return_internal_reference<>())
    .def("GetUserEventAction",    &G4RunManager::GetUserEventAction,
   return_internal_reference<>())
    .def("GetUserStackingAction", &G4RunManager::GetUserStackingAction,
   return_internal_reference<>())
    .def("GetUserTrackingAction", &G4RunManager::GetUserTrackingAction,
   return_internal_reference<>())
    .def("GetUserSteppingAction", &G4RunManager::GetUserSteppingAction,
   return_internal_reference<>())
    // ---
    .def("AbortRun",             &G4RunManager::AbortRun,
   f_AbortRun((arg("soft_abort")=false),
        "Abort run (event loop)."))
    .def("AbortEvent",           &G4RunManager::AbortEvent)
    .def("DefineWorldVolume",    &G4RunManager::DefineWorldVolume,
                                 f_DefineWorldVolume())
    .def("DumpRegion",           f1_DumpRegion)
    .def("DumpRegion",           f2_DumpRegion, f_DumpRegion())
    .def("rndmSaveThisRun",      &G4RunManager::rndmSaveThisRun)
    .def("rndmSaveThisEvent",    &G4RunManager::rndmSaveThisEvent)
    .def("RestoreRandomNumberStatus",
                                 &G4RunManager::RestoreRandomNumberStatus)
    .def("SetRandomNumberStore", &G4RunManager::SetRandomNumberStore)
    .def("GetRandomNumberStore", &G4RunManager::GetRandomNumberStore)
    .def("SetRandomNumberStoreDir", &G4RunManager::SetRandomNumberStoreDir)
    .def("GeometryHasBeenModified", &G4RunManager::GeometryHasBeenModified)
    .def("PhysicsHasBeenModified",  &G4RunManager::PhysicsHasBeenModified)
    .def("GetGeometryToBeOptimized",&G4RunManager::GetGeometryToBeOptimized)
    .def("GetCurrentRun",  &G4RunManager::GetCurrentRun,
    return_value_policy<reference_existing_object>())
    .def("GetCurrentEvent", &G4RunManager::GetCurrentEvent,
    return_value_policy<reference_existing_object>())
    .def("SetRunIDCounter",        &G4RunManager::SetRunIDCounter)

#if G4VERSION_NUMBER >= 932
    .def("GetVersionString",     &G4RunManager::GetVersionString,
    return_value_policy<reference_existing_object>())
    .def("GetRandomNumberStoreDir", &G4RunManager::GetRandomNumberStoreDir,
    return_internal_reference<>())
#else
    .def("GetVersionString",        &G4RunManager::GetVersionString)
    .def("GetRandomNumberStoreDir", &G4RunManager::GetRandomNumberStoreDir)
#endif
    ;

    // reduced functionality...
    // void SetPrimaryTransformer(G4PrimaryTransformer* pt)
    // void SetNumberOfAdditionalWaitingStacks(G4int iAdd)
    // void CutOffHasBeenModified()
    // void SetGeometryToBeOptimized(G4bool vl)
    // const G4Event* GetPreviousEvent(G4int i) const
    // void SetNumberOfEventsToBeStored(G4int val)
    // void SetDCtable(G4DCtable* DCtbl)

}
