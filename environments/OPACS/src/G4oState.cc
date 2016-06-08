// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4oState.cc,v 1.6 2000/03/17 08:57:48 barrand Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//#define DEBUG

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
G4oState::~G4oState(
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
G4bool G4oState::Notify (
 G4ApplicationState requestedState
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4StateManager* statM = G4StateManager::GetStateManager();
  G4ApplicationState previousState = statM->GetPreviousState();

  if(previousState==Idle && requestedState==GeomClosed) {             
    //beginOfRun
    G4String string   = name; 
    string += "RunBegin.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
  } else if(previousState==GeomClosed && requestedState==Idle) {      
    //endOfRun
    G4String string = name; 
    string += "RunEnd.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
  } else if(previousState==GeomClosed && requestedState==EventProc) { 
    //beginOfEvent
    G4String string = name; 
    string += "EventBegin.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
  } else if(previousState==EventProc && requestedState==GeomClosed) { 
    // EndOfEvent
    G4String string = name; 
    string += "EventEnd.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
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









