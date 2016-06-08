// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4oState.cc,v 1.2 1999/04/16 10:03:39 barrand Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
//#define DEBUG

//#include "G4ios.hh"

//G4
#include <globals.hh>
#include <G4UImanager.hh>
#include <G4StateManager.hh>

//OPACS
#include <CPrinter.h>
#include <Wo.h>

#include <G4o.h>
#include <G4oState.hh>

/***************************************************************************/
G4oState::G4oState(
 G4String a_name
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  name = a_name;
  G4UImanager::GetUIpointer()->SetSession(this);  //So that Pause works..
}
/***************************************************************************/
G4bool G4oState::Notify (
 G4ApplicationState requestedState
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4StateManager*            statM = G4StateManager::GetStateManager();
  G4ApplicationState previousState = statM->GetPreviousState();

  if(previousState==Idle && requestedState==GeomClosed) {             //beginOfRun
    G4String string   = "RunBegin.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)(name + string).data());
  } else if(previousState==GeomClosed && requestedState==Idle) {      //endOfRun
    G4String string   = "RunEnd.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)(name + string).data());
  } else if(previousState==GeomClosed && requestedState==EventProc) { //beginOfEvent
    G4String string   = "EventBegin.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)(name + string).data());
  } else if(previousState==EventProc && requestedState==GeomClosed) { // EndOfEvent
    G4String string   = "EventEnd.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)(name + string).data());
  }

  return true;
}
/***************************************************************************/
void G4oState::PauseSessionStart (
 G4String a_message
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_message=="G4_pause> ") { 
    WoProcessEvents ();
  }
}









