// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StackManager.cc,v 1.2 1999-12-15 14:49:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 09/Dec/96 M.Asai
//

#include "G4StackManager.hh"
#include "G4StackingMessenger.hh"
#include "evmandefs.hh"
#include "G4ios.hh"

G4StackManager::G4StackManager()
:verboseLevel(0),userStackingAction(NULL)
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
}

const G4StackManager & G4StackManager::operator=
(const G4StackManager &) { return *this; }
int G4StackManager::operator==(const G4StackManager &) 
const{ return false; }
int G4StackManager::operator!=(const G4StackManager &) 
const{ return true; }

G4int G4StackManager::PushOneTrack(G4Track *newTrack)
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
  }
  else
  {
    G4StackedTrack * newStackedTrack = new G4StackedTrack( newTrack );
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
    }
  }

  return GetNUrgentTrack();
}


G4Track * G4StackManager::PopNextTrack()
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
    if(userStackingAction) userStackingAction->NewStage();
#ifdef G4VERBOSE
    if( verboseLevel > 0 ) G4cout << "     " << GetNUrgentTrack()
                      << " urgent tracks and " << GetNWaitingTrack()
                      << " waiting tracks." << G4endl;
#endif
    if( ( GetNUrgentTrack()==0 ) && ( GetNWaitingTrack()==0 ) ) return NULL;
  }

  G4StackedTrack * selectedStackedTrack = urgentStack->PopFromStack();
  G4Track * selectedTrack = selectedStackedTrack->GetTrack();

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
  while( (aStackedTrack=tmpStack.PopFromStack()) != NULL )
  {
    G4ClassificationOfNewTrack classification = 
      userStackingAction->ClassifyNewTrack( aStackedTrack->GetTrack() );
    switch (classification)
    {
      case fKill:
        delete aStackedTrack->GetTrack();
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

    while( (aStackedTrack=tmpStack.PopFromStack()) != NULL )
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
        }
      }
    }
  }

  return n_passedFromPrevious;
}


