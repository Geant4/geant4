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
// $Id: G4StackManager.cc,v 1.6 2001-07-19 00:14:17 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 09/Dec/96 M.Asai
//

#include "G4StackManager.hh"
#include "G4StackingMessenger.hh"
#include "G4VTrajectory.hh"
#include "evmandefs.hh"
#include "G4ios.hh"

G4StackManager::G4StackManager()
:userStackingAction(0),verboseLevel(0),numberOfAdditionalWaitingStacks(0)
{
  theMessenger = new G4StackingMessenger(this);
  urgentStack = new G4TrackStack;
  waitingStack = new G4TrackStack;
  postponeStack = new G4TrackStack;
}

G4StackManager::~G4StackManager()
{
  if(userStackingAction) delete userStackingAction;
  delete urgentStack;
  delete waitingStack;
  delete postponeStack;
  delete theMessenger;
  if(numberOfAdditionalWaitingStacks>0) {
    for(int i=0;i<numberOfAdditionalWaitingStacks;i++) {
      delete additionalWaitingStacks[i];
    }
  }
}

const G4StackManager & G4StackManager::operator=
(const G4StackManager &) { return *this; }
G4int G4StackManager::operator==(const G4StackManager &) 
const{ return false; }
G4int G4StackManager::operator!=(const G4StackManager &) 
const{ return true; }

G4int G4StackManager::PushOneTrack(G4Track *newTrack,G4VTrajectory *newTrajectory)
{
  G4ClassificationOfNewTrack classification;
  if(userStackingAction) 
  { classification = userStackingAction->ClassifyNewTrack( newTrack ); }
  else
  { classification = DefaultClassification( newTrack ); }

  if(classification==fKill)   // delete newTrack without stacking
  {
#ifdef G4VERBOSE
    if( verboseLevel > 0 )
    {
      G4cout << "   ---> G4Track " << newTrack << " (trackID "
	 << newTrack->GetTrackID() << ", parentID "
	 << newTrack->GetParentID() << ") is not to be stored." << G4endl;
    }
#endif
    delete newTrack;
    delete newTrajectory;
  }
  else
  {
    G4StackedTrack * newStackedTrack = new G4StackedTrack( newTrack, newTrajectory );
    switch (classification)
    {
      case fUrgent:
        urgentStack->PushToStack( newStackedTrack );
        break;
      case fWaiting:
        waitingStack->PushToStack( newStackedTrack );
        break;
      case fPostpone:
        postponeStack->PushToStack( newStackedTrack );
        break;
      case fKill:
        break;
      default:
        G4int i = classification - 10;
        if(i<1||i>numberOfAdditionalWaitingStacks) {
          G4Exception("G4StackManager : invalid classification");
        } else {
          additionalWaitingStacks[i-1]->PushToStack( newStackedTrack );
        }
        break;
    }
  }

  return GetNUrgentTrack();
}


G4Track * G4StackManager::PopNextTrack(G4VTrajectory**newTrajectory)
{
#ifdef G4VERBOSE
  if( verboseLevel > 0 )
  {
    G4cout << "### pop requested out of " 
         << GetNUrgentTrack() << " stacked tracks." << G4endl;
  }
#endif

  while( GetNUrgentTrack() == 0 )
  {
#ifdef G4VERBOSE
    if( verboseLevel > 0 ) G4cout << "### " << GetNWaitingTrack()
                      << " waiting tracks are re-classified to" << G4endl;
#endif
    waitingStack->TransferTo(urgentStack);
    if(numberOfAdditionalWaitingStacks>0) {
      for(int i=0;i<numberOfAdditionalWaitingStacks;i++) {
        if(i==0) {
          additionalWaitingStacks[0]->TransferTo(waitingStack);
        } else {
          additionalWaitingStacks[i]->TransferTo(additionalWaitingStacks[i-1]);
        }
      }
    }
    if(userStackingAction) userStackingAction->NewStage();
#ifdef G4VERBOSE
    if( verboseLevel > 0 ) G4cout << "     " << GetNUrgentTrack()
                      << " urgent tracks and " << GetNWaitingTrack()
                      << " waiting tracks." << G4endl;
#endif
    if( ( GetNUrgentTrack()==0 ) && ( GetNWaitingTrack()==0 ) ) return 0;
  }

  G4StackedTrack * selectedStackedTrack = urgentStack->PopFromStack();
  G4Track * selectedTrack = selectedStackedTrack->GetTrack();
  *newTrajectory = selectedStackedTrack->GetTrajectory();

#ifdef G4VERBOSE
  if( verboseLevel > 1 )
  {
    G4cout << "Selected G4StackedTrack : " << selectedStackedTrack 
         << " with G4Track " << selectedStackedTrack->GetTrack()
	 << " (trackID " << selectedStackedTrack->GetTrack()->GetTrackID()
	 << ", parentID " << selectedStackedTrack->GetTrack()->GetParentID()
	 << ")" << G4endl;
  }
#endif

  delete selectedStackedTrack;
  return selectedTrack;
}

void G4StackManager::ReClassify()
{
  G4StackedTrack * aStackedTrack;
  G4TrackStack tmpStack;

  if( !userStackingAction ) return;
  if( GetNUrgentTrack() == 0 ) return;

  urgentStack->TransferTo(&tmpStack);
  while( (aStackedTrack=tmpStack.PopFromStack()) != 0 )
  {
    G4ClassificationOfNewTrack classification = 
      userStackingAction->ClassifyNewTrack( aStackedTrack->GetTrack() );
    switch (classification)
    {
      case fKill:
        delete aStackedTrack->GetTrack();
        delete aStackedTrack->GetTrajectory();
        delete aStackedTrack;
        break;
      case fUrgent:
        urgentStack->PushToStack( aStackedTrack );
        break;
      case fWaiting:
        waitingStack->PushToStack( aStackedTrack );
        break;
      case fPostpone:
        postponeStack->PushToStack( aStackedTrack );
        break;
      default:
        G4int i = classification - 10;
        if(i<1||i>numberOfAdditionalWaitingStacks) {
          G4Exception("G4StackManager : invalid classification");
        } else {
          additionalWaitingStacks[i-1]->PushToStack( aStackedTrack );
        }
        break;
    }
  }
}

G4int G4StackManager::PrepareNewEvent()
{
  if(userStackingAction) userStackingAction->PrepareNewEvent();

  G4int n_passedFromPrevious = 0;

  if( GetNPostponedTrack() > 0 )
  {
#ifdef G4VERBOSE
    if( verboseLevel > 0 )
    {
      G4cout << GetNPostponedTrack() 
           << " postponed tracked are now shifted to the stack." << G4endl;
    }
#endif

    G4StackedTrack * aStackedTrack;
    G4TrackStack tmpStack;

    postponeStack->TransferTo(&tmpStack);

    while( (aStackedTrack=tmpStack.PopFromStack()) != 0 )
    {
      G4Track* aTrack = aStackedTrack->GetTrack();
      aTrack->SetParentID(-1);
      G4ClassificationOfNewTrack classification;
      if(userStackingAction) 
      { classification = userStackingAction->ClassifyNewTrack( aTrack ); }
      else
      { classification = DefaultClassification( aTrack ); }

      if(classification==fKill)
      {
        delete aTrack;
        delete aStackedTrack;
      }
      else
      {
        aTrack->SetTrackID(++n_passedFromPrevious);
        switch (classification)
        {
          case fUrgent:
            urgentStack->PushToStack( aStackedTrack );
            break;
          case fWaiting:
            waitingStack->PushToStack( aStackedTrack );
            break;
          case fPostpone:
            postponeStack->PushToStack( aStackedTrack );
            break;
          case fKill:
            break;
          default:
            G4int i = classification - 10;
            if(i<1||i>numberOfAdditionalWaitingStacks) {
              G4Exception("G4StackManager : invalid classification");
            } else {
              additionalWaitingStacks[i-1]->PushToStack( aStackedTrack );
            }
            break;
        }
      }
    }
  }

  return n_passedFromPrevious;
}

void G4StackManager::SetNumberOfAdditionalWaitingStacks(G4int iAdd)
{
  if(iAdd > numberOfAdditionalWaitingStacks)
  {
    for(int i=numberOfAdditionalWaitingStacks;i<iAdd;i++)
    {
      G4TrackStack* newStack = new G4TrackStack;
      additionalWaitingStacks.push_back(newStack);
    }
    numberOfAdditionalWaitingStacks = iAdd;
  }
  else if (iAdd < numberOfAdditionalWaitingStacks)
  {
    for(int i=numberOfAdditionalWaitingStacks;i>iAdd;i--)
    {
      delete additionalWaitingStacks[i];
    }
  }
}

void G4StackManager::TransferStackedTracks(G4ClassificationOfNewTrack origin, G4ClassificationOfNewTrack destination)
{
  G4TrackStack* originStack = 0;
  switch(origin)
  {
    case fUrgent:
      originStack = urgentStack;
      break;
    case fWaiting:
      originStack = waitingStack;
      break;
    case fPostpone:
      originStack = postponeStack;
      break;
    case fKill:
      break;
    default:
      int i = origin - 10;
      if(i<=numberOfAdditionalWaitingStacks) originStack = additionalWaitingStacks[i-1];
      break;
  }
  if(!originStack) return;

  if(destination==fKill)
  {
    G4StackedTrack * aStackedTrack;
    while( (aStackedTrack=originStack->PopFromStack()) != 0 )
    {
      delete aStackedTrack->GetTrack();
      delete aStackedTrack;
    }
  } 
  else
  {
    G4TrackStack* targetStack = 0;
    switch(destination)
    {
      case fUrgent:
        targetStack = urgentStack;
        break;
      case fWaiting:
        targetStack = waitingStack;
        break;
      case fPostpone:
        targetStack = postponeStack;
        break;
      default:
        int i = origin - 10;
        if(i<=numberOfAdditionalWaitingStacks) targetStack = additionalWaitingStacks[i-1];
        break;
    }
    if(!targetStack) return;
    if(originStack==targetStack) return;
    originStack->TransferTo(targetStack);
  }
  return;
}





