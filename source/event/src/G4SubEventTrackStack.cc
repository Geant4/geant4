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
// G4SubEventTrackStack class implementation
//
// Author: Makoto Asai (JLab) - 23/Aug/23
// --------------------------------------------------------------------

#include "G4SubEventTrackStack.hh"
#include "G4StackedTrack.hh"
#include "G4Track.hh"
#include "G4VTrajectory.hh"
#include "G4Event.hh"

G4SubEventTrackStack::~G4SubEventTrackStack()
{
  clearAndDestroy();
}

void G4SubEventTrackStack::clearAndDestroy()
{
  if(fCurrentSE!=nullptr) 
  {
    fCurrentSE->clearAndDestroy();
    delete fCurrentSE;
    fCurrentSE = nullptr;
  }
}

void G4SubEventTrackStack::PrepareNewEvent(G4Event* ev)
{
  if(fCurrentSE!=nullptr)
  {
    G4ExceptionDescription ed;
    ed << fCurrentSE->size() << " sub-events still remains in the previous event. PANIC!!!";
    G4Exception("G4SubEventTrackStack::PrepareNewEvent()","SubEvt7001",FatalException,ed);
  }
  fCurrentSE = nullptr;
  fCurrentEvent = ev;
}
    
void G4SubEventTrackStack::PushToStack(const G4StackedTrack& aStackedTrack)
{
  if(fCurrentSE==nullptr)
  {
    fCurrentSE = new G4SubEvent(fSubEventType,fMaxEnt); 
  }
  else if(fCurrentSE->size()==fMaxEnt)
  {
    // current sus-event is already full. Transfer it to G4Event and create a new sub-event.
    auto nSubEv = fCurrentEvent->StoreSubEvent(fSubEventType,fCurrentSE);
    if(verboseLevel>1)
    {
      G4cout << "### event id " << fCurrentEvent->GetEventID()
             << " -- sub-evnet " << nSubEv << " with " << fCurrentSE->size()
             << " tracks is stored" << G4endl;
    }
    fCurrentSE = new G4SubEvent(fSubEventType,fMaxEnt);
  }
  fCurrentSE->PushToStack(aStackedTrack);
}

void G4SubEventTrackStack::ReleaseSubEvent()
{
  // This method should be invoked at the end of processing an event.
  // A sub-event that is not yet full is transferred to G4Event for
  // post-event-loop processing.
  if(fCurrentEvent==nullptr)
  {
    G4Exception("G4SubEventTrackStack::ReleaseSubEvent()","SubEvt7002",FatalException,
                "Valid pointer of the current event is not set. PANIC!!");
    return; // NOLINT: Explicit return required to silence Coverity
  }
  if(fCurrentSE!=nullptr)
  {
    auto nSubEv = fCurrentEvent->StoreSubEvent(fSubEventType,fCurrentSE);
    if(verboseLevel>1)
    {
      G4cout << "### event id " << fCurrentEvent->GetEventID()
             << " -- sub-evnet " << nSubEv << " with " << fCurrentSE->size()
             << " tracks is stored" << G4endl;
    }
    fCurrentSE = nullptr;
  }
  fCurrentEvent = nullptr;
}

//G4StackedTrack G4SubEventTrackStack::PopFromStack()
//{
//  G4Exception("G4SubEventTrackStack::PopFromStack","SubEvt7000",
//       FatalException,"This method must not be invoked.");
//  return G4StackedTrack(nullptr);
//}

