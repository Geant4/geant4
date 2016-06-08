// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4oState.cc,v 1.3.2.1.2.1 1999/12/07 20:46:52 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
//#define DEBUG

// GB : put this include before the GL ones. 
// It avoid a clash with STL includes on Linux.
#include "g4rw/tvhdict.h"

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









