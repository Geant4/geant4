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
// $Id: ParRunManager.cc,v 1.1 2002-03-05 15:22:18 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------

#ifdef G4USE_TOPC

// G4Timer.hh has to be loaded *before* globals.hh...
//
#include "G4Timer.hh"

#include "G4RunManager.hh"

#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4VPersistencyManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"

#include "ParRunManager.hh"
#include "ParMarshaledObj.hh"
#include "ParRandomState.hh"

#include "G4ios.hh"
#include "g4std/strstream"

// Define static data members of ParRunManager
// They are set by ImportDoEventLoopLocals
G4StateManager* ParRunManager::stateManager;
G4int ParRunManager::n_event;
G4int ParRunManager::n_select;
G4String ParRunManager::msg;
ParRunManager* ParRunManager::myRunManager = 0;

// TOP-C callbacks (static functions)
//   If this->*(&ParRunManager::GenerateEventInput), created a thunk,
//   it would be ideal here.
TOPC_BUF ParRunManager::MyGenerateEventInput()
  { return myRunManager->GenerateEventInput(); }
TOPC_BUF ParRunManager::MyDoEvent( void *input_buf )
  { return myRunManager->DoEvent(input_buf); }
TOPC_ACTION ParRunManager::MyCheckEventResult( void * input_buf, void *buf )
  { return myRunManager->CheckEventResult(input_buf, buf); }

static void trace_event_input( void *input )
{ G4cout << "Event " << *(G4int *)input << G4endl; }

void ParRunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  TOPC_OPT_trace_input = trace_event_input;

  G4StateManager* stateManager = G4StateManager::GetStateManager();

#ifdef G4_ORIGINAL
  cout << "ParRunManager::DoEventLoop" << endl;
  G4RunManager::DoEventLoop(n_event, macroFile, n_select);
  return;
#endif

  if(verboseLevel>0) 
  { timer->Start(); }

  G4String msg;
  if(macroFile!=0)
  { 
    if(n_select<0) n_select = n_event;
    msg = "/control/execute ";
    msg += macroFile;
  }
  else
  { n_select = -1; }

  // BeginOfEventAction() and EndOfEventAction() would normally be
  // called inside G4EventManager::ProcessOneEvent() in ParRunManager::DoEvent
  // on slave.  Since this is often where hits are collected and where
  // histogram data is first stored, must do it
  // on master instead, to process hits.  So, set user event action to NULL.
  // Keep private copy, and execute it in CheckTaskResult().
  // If user later does:  SetUserAction( (G4UserEventAction *)NULL );
  //   we won't notice and copy it to origUserEventAction.
  // How paranoid do we want to be?  Should we create trvialUserEventAction?
  if ( eventManager->GetUserEventAction() )
  {
    origUserEventAction = eventManager->GetUserEventAction();
    SetUserAction( (G4UserEventAction *)0 );
  }

  // Make these variables accessible to TOP-C callback functions
  ImportDoEventLoopLocals( stateManager, n_event, n_select, msg );

  // This is where all the parallelism occurs
  TOPC_master_slave(MyGenerateEventInput, MyDoEvent, MyCheckEventResult, 0);

  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
/*
    { G4cout << "  Run Aborted after " << i_event << " events processed." << G4endl; }
*/
    { G4cout << "  Run Aborted." << G4endl; }
    else
    { G4cout << "  Number of events processed : " << n_event << G4endl; }
    G4cout << "  "  << *timer << G4endl;
  }

}

TOPC_BUF ParRunManager::GenerateEventInput()
{
  static G4int i_event = -1;

  if(runAborted) return NOTASK;

  i_event++;
  if (i_event >= n_event) return NOTASK;
  else
  { ParMarshaledRandomState buf = ParMarshaledRandomState( i_event );
    return TOPC_MSG( buf.GetBuffer(), buf.BufferSize() );
  }
  // else return TOPC_MSG( &i_event, sizeof(i_event) );
}

TOPC_BUF ParRunManager::DoEvent( void *input_buf )
{
  G4int i_event;
  ParMarshaledRandomState buf = ParMarshaledRandomState( input_buf );
  buf.unmarshalEventIDandSetState( i_event );

    stateManager->SetNewState(EventProc);

    currentEvent = GenerateEvent(i_event);

    eventManager->ProcessOneEvent(currentEvent);

    MarshaledObj *aMarshaledObj = new MarshaledHCofThisEvent(currentEvent);
    return TOPC_MSG(aMarshaledObj->GetBuffer(), aMarshaledObj->BufferSize());
}

TOPC_ACTION ParRunManager::CheckEventResult( void * input_buf, void *output_buf )
{
  G4int i_event;
  MarshaledObj buf = MarshaledObj(input_buf);
  buf.Unmarshal(i_event);

  stateManager->SetNewState(EventProc);
  if(!userPrimaryGeneratorAction)
  {
    G4Exception
    ("G4RunManager::BeamOn - G4VUserPrimaryGeneratorAction is not defined.");
  }

  //This creates a trivial event in lieu of GenerateEvent(i_event);
  currentEvent = new G4Event(i_event);

  //Original UserEventAction was saved and set to NULL.  Do it now on master.
  if (origUserEventAction)
    origUserEventAction->BeginOfEventAction( currentEvent );

  // Patches minor bug in Geant4 4.0;  see docs/BUG.userStackingAction
  // When Geant4 4.0 sees empty event, it still calls userStackingAction.
  // On master, only trivial events exist, so we delete userStackingAction
  SetUserAction( (G4UserStackingAction*)0 );
  eventManager->ProcessOneEvent(currentEvent); // Processing the trivial event

  // Called with output_buf and no size, creates object for unmarshaling
  MarshaledHCofThisEvent marshaledObj(output_buf);
  marshaledObj.UnmarshalSlaveHCofThisEvent();
  // UnmarshalSlaveHCofThisEvent placed slave Hits in ParRunManager::currentEvent

  // Original UserEventAction was saved and set to NULL.  Do it now on master.
  if (origUserEventAction)
    origUserEventAction->EndOfEventAction( currentEvent );

  AnalyzeEvent(currentEvent);

    if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msg);
    stateManager->SetNewState(GeomClosed);
    StackPreviousEvent(currentEvent);
    currentEvent = 0;
    //Move this to GenerateEventInput:
    // if(runAborted) break;
    return NO_ACTION;
}

#endif
