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
//
// $Id: G4oState.cc,v 1.8 2002/12/06 16:19:22 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-02 $
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

  if(previousState==G4State_Idle && requestedState==G4State_GeomClosed) {             
    //beginOfRun
    G4String string   = name; 
    string += "RunBegin.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
  } else if(previousState==G4State_GeomClosed && requestedState==G4State_Idle) {      
    //endOfRun
    G4String string = name; 
    string += "RunEnd.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
  } else if(previousState==G4State_GeomClosed && requestedState==G4State_EventProc) { 
    //beginOfEvent
    G4String string = name; 
    string += "EventBegin.osh"; 
    OShellExecuteFile (G4oGetShell(),(char*)string.data());
  } else if(previousState==G4State_EventProc && requestedState==G4State_GeomClosed) { 
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









