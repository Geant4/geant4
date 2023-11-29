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
// G4StackManager class implementation
//
// Author: Makoto Asai, 1996
// --------------------------------------------------------------------

#include "G4StackManager.hh"
#include "G4StackingMessenger.hh"
#include "G4VTrajectory.hh"
#include "G4ios.hh"

#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

// Needed for temporal service
//
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

G4StackManager::G4StackManager()
{
  theMessenger = new G4StackingMessenger(this);
#ifdef G4_USESMARTSTACK
  urgentStack = new G4SmartTrackStack;
  // G4cout << "+++ G4StackManager uses G4SmartTrackStack. +++" << G4endl;
#else
  urgentStack = new G4TrackStack(5000);
  // G4cout << "+++ G4StackManager uses ordinary G4TrackStack. +++" << G4endl;
#endif
  waitingStack = new G4TrackStack(1000);
  postponeStack = new G4TrackStack(1000);
}

G4StackManager::~G4StackManager()
{
  delete userStackingAction; 

#ifdef G4VERBOSE
  if(verboseLevel>0)
  {
    G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << G4endl;
    G4cout << " Maximum number of tracks in the urgent stack : " << urgentStack->GetMaxNTrack() << G4endl;
    G4cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << G4endl;
  }
#endif
  delete urgentStack;
  delete waitingStack;
  delete postponeStack;
  delete theMessenger;
  if(numberOfAdditionalWaitingStacks>0)
  {
    for(G4int i=0; i<numberOfAdditionalWaitingStacks; ++i)
    {
      delete additionalWaitingStacks[i];
    }
  }
}

G4int G4StackManager::
PushOneTrack(G4Track* newTrack, G4VTrajectory* newTrajectory)
{
  const G4ParticleDefinition* pd = newTrack->GetParticleDefinition();
  if(pd->GetParticleDefinitionID() < 0)
  {
    G4ExceptionDescription ED;
    ED << "A track without proper process manager is pushed \
           into the track stack.\n"
       << " Particle name : " << pd->GetParticleName() << " -- ";
    if(newTrack->GetParentID()<0)
    {
      ED << "created by a primary particle generator.";
    }
    else
    { 
      const G4VProcess* vp = newTrack->GetCreatorProcess();
      if(vp != nullptr)
      {
        ED << "created by " << vp->GetProcessName() << ".";
      }
      else
      {
        ED << "creaded by unknown process.";
      }
    }
    G4Exception("G4StackManager::PushOneTrack","Event10051",
                 FatalException,ED);
    delete newTrack;
    return GetNUrgentTrack();
  }
    
  G4ClassificationOfNewTrack classification = DefaultClassification( newTrack );
  if(userStackingAction != nullptr) 
  {
    classification = userStackingAction->ClassifyNewTrack( newTrack );
  }

  if(classification==fKill)   // delete newTrack without stacking
  {
#ifdef G4VERBOSE
    if( verboseLevel > 1 )
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
    G4StackedTrack newStackedTrack( newTrack, newTrajectory );
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
      default:
        G4int i = classification - 10;
        if(i<1 || i>numberOfAdditionalWaitingStacks)
        {
          G4ExceptionDescription ED;
          ED << "invalid classification " << classification << G4endl;
          G4Exception("G4StackManager::PushOneTrack", "Event0051",
                      FatalException,ED);
        }
        else
        {
          additionalWaitingStacks[i-1]->PushToStack( newStackedTrack );
        }
        break;
    }
  }
  return GetNUrgentTrack();
}


G4Track* G4StackManager::PopNextTrack(G4VTrajectory** newTrajectory)
{
#ifdef G4VERBOSE
  if( verboseLevel > 1 )
  {
    G4cout << "### pop requested out of " 
           << GetNUrgentTrack() << " stacked tracks." << G4endl;
  }
#endif

  while( GetNUrgentTrack() == 0 )
  {
#ifdef G4VERBOSE
    if( verboseLevel > 1 )
    {
      G4cout << "### " << GetNWaitingTrack()
             << " waiting tracks are re-classified to" << G4endl;
    }
#endif
    waitingStack->TransferTo(urgentStack);
    if(numberOfAdditionalWaitingStacks>0)
    {
      for(G4int i=0; i<numberOfAdditionalWaitingStacks; ++i)
      {
        if(i==0)
        {
          additionalWaitingStacks[0]->TransferTo(waitingStack);
        }
        else
        {
          additionalWaitingStacks[i]->TransferTo(additionalWaitingStacks[i-1]);
        }
      }
    }
    if(userStackingAction != nullptr)
    {
      userStackingAction->NewStage();
    }

#ifdef G4VERBOSE
    if( verboseLevel > 1 )
      G4cout << "     " << GetNUrgentTrack()
             << " urgent tracks and " << GetNWaitingTrack()
             << " waiting tracks." << G4endl;
#endif
    if( ( GetNUrgentTrack()==0 ) && ( GetNWaitingTrack()==0 ) )
      return nullptr;
  }

  G4StackedTrack selectedStackedTrack = urgentStack->PopFromStack();
  G4Track * selectedTrack = selectedStackedTrack.GetTrack();
  *newTrajectory = selectedStackedTrack.GetTrajectory();

#ifdef G4VERBOSE
  if( verboseLevel > 2 )
  {
    G4cout << "Selected G4StackedTrack : " << &selectedStackedTrack
           << " with G4Track " << selectedStackedTrack.GetTrack()
           << " (trackID " << selectedStackedTrack.GetTrack()->GetTrackID()
           << ", parentID " << selectedStackedTrack.GetTrack()->GetParentID()
           << ")" << G4endl;
  }
#endif

  return selectedTrack;
}

void G4StackManager::ReClassify()
{
  G4StackedTrack aStackedTrack;
  G4TrackStack tmpStack;
  
  if( userStackingAction == nullptr ) return;
  if( GetNUrgentTrack() == 0 ) return;
  
  urgentStack->TransferTo(&tmpStack);
  while( tmpStack.GetNTrack() > 0 )
  {
    aStackedTrack=tmpStack.PopFromStack();
    G4ClassificationOfNewTrack classification =
    userStackingAction->ClassifyNewTrack( aStackedTrack.GetTrack() );
    switch (classification)
    {
      case fKill:
        delete aStackedTrack.GetTrack();
        delete aStackedTrack.GetTrajectory();
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
        if(i<1||i>numberOfAdditionalWaitingStacks)
        {
          G4ExceptionDescription ED;
          ED << "invalid classification " << classification << G4endl;
          G4Exception("G4StackManager::ReClassify", "Event0052",
                      FatalException, ED);
        }
        else
        {
          additionalWaitingStacks[i-1]->PushToStack( aStackedTrack );
        }
        break;
    }
  }
}

G4int G4StackManager::PrepareNewEvent()
{
  if(userStackingAction != nullptr)
  {
    userStackingAction->PrepareNewEvent();
  }
  
  // Set the urgentStack in a defined state. Not doing it would
  // affect reproducibility
  //
  urgentStack->clearAndDestroy();
  
  G4int n_passedFromPrevious = 0;
  
  if( GetNPostponedTrack() > 0 )
  {
#ifdef G4VERBOSE
    if( verboseLevel > 1 )
    {
      G4cout << GetNPostponedTrack()
             << " postponed tracked are now shifted to the stack." << G4endl;
    }
#endif
    
    G4StackedTrack aStackedTrack;
    G4TrackStack   tmpStack;
    
    postponeStack->TransferTo(&tmpStack);
    
    while( tmpStack.GetNTrack() > 0 )
    {
      aStackedTrack=tmpStack.PopFromStack();
      G4Track* aTrack = aStackedTrack.GetTrack();
      aTrack->SetParentID(-1);
      G4ClassificationOfNewTrack classification;
      if(userStackingAction != nullptr)
      {
        classification = userStackingAction->ClassifyNewTrack( aTrack );
      }
      else
      {
        classification = DefaultClassification( aTrack );
      }
      
      if(classification==fKill)
      {
        delete aTrack;
        delete aStackedTrack.GetTrajectory();
      }
      else
      {
        aTrack->SetTrackID(-(++n_passedFromPrevious));
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
          default:
            G4int i = classification - 10;
            if(i<1||i>numberOfAdditionalWaitingStacks)
            {
              G4ExceptionDescription ED;
              ED << "invalid classification " << classification << G4endl;
              G4Exception("G4StackManager::PrepareNewEvent", "Event0053",
                          FatalException, ED);
            }
            else
            {
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
    for(G4int i=numberOfAdditionalWaitingStacks; i<iAdd; ++i)
    {
      auto* newStack = new G4TrackStack;
      additionalWaitingStacks.push_back(newStack);
    }
    numberOfAdditionalWaitingStacks = iAdd;
  }
  else if (iAdd < numberOfAdditionalWaitingStacks)
  {
    for(G4int i=numberOfAdditionalWaitingStacks; i>iAdd; --i)
    {
      delete additionalWaitingStacks[i];
    }
  }
}

void G4StackManager::
TransferStackedTracks(G4ClassificationOfNewTrack origin,
                      G4ClassificationOfNewTrack destination)
{
  if(origin==destination) return;
  if(origin==fKill) return;
  G4TrackStack* originStack = nullptr;
  switch(origin)
  {
    case fUrgent:
      originStack = nullptr;
      break;
    case fWaiting:
      originStack = waitingStack;
      break;
    case fPostpone:
      originStack = postponeStack;
      break;
    default:
      G4int i = origin - 10;
      if(i<=numberOfAdditionalWaitingStacks)
      {
        originStack = additionalWaitingStacks[i-1];
      }
      break;
  }
  
  if(destination==fKill)
  {
    if(originStack != nullptr)
    {
      originStack->clearAndDestroy();
    }
    else
    {
      urgentStack->clearAndDestroy();
    }
  }
  else
  {
    G4TrackStack* targetStack = nullptr;
    switch(destination)
    {
      case fUrgent:
        targetStack = nullptr;
        break;
      case fWaiting:
        targetStack = waitingStack;
        break;
      case fPostpone:
        targetStack = postponeStack;
        break;
      default:
        G4int i = destination - 10;
        if(i<=numberOfAdditionalWaitingStacks)
        {
          targetStack = additionalWaitingStacks[i-1];
        }
        break;
    }
    if(originStack != nullptr)
    {
      if(targetStack != nullptr)
      {
        originStack->TransferTo(targetStack);
      }
      else
      {
        originStack->TransferTo(urgentStack);
      }
    }
    else
    {
      urgentStack->TransferTo(targetStack);
    }
  }
  return;
}

void G4StackManager::
TransferOneStackedTrack(G4ClassificationOfNewTrack origin,
                        G4ClassificationOfNewTrack destination)
{
  if(origin==destination) return;
  if(origin==fKill) return;
  G4TrackStack* originStack = nullptr;
  switch(origin)
  {
    case fUrgent:
      originStack = nullptr;
      break;
    case fWaiting:
      originStack = waitingStack;
      break;
    case fPostpone:
      originStack = postponeStack;
      break;
    default:
      G4int i = origin - 10;
      if(i<=numberOfAdditionalWaitingStacks)
      {
        originStack = additionalWaitingStacks[i-1];
      }
      break;
  }
  
  G4StackedTrack aStackedTrack;
  if(destination==fKill)
  {
    if( originStack != nullptr && (originStack->GetNTrack() != 0u) )
    {
      aStackedTrack = originStack->PopFromStack();
      delete aStackedTrack.GetTrack();
      delete aStackedTrack.GetTrajectory();
    }
    else if (urgentStack->GetNTrack() != 0u )
    {
      aStackedTrack = urgentStack->PopFromStack();
      delete aStackedTrack.GetTrack();
      delete aStackedTrack.GetTrajectory();
    }
  }
  else
  {
    G4TrackStack* targetStack = nullptr;
    switch(destination)
    {
      case fUrgent:
        targetStack = nullptr;
        break;
      case fWaiting:
        targetStack = waitingStack;
        break;
      case fPostpone:
        targetStack = postponeStack;
        break;
      default:
        G4int i = destination - 10;
        if(i<=numberOfAdditionalWaitingStacks)
        {
          targetStack = additionalWaitingStacks[i-1];
        }
        break;
    }
    if((originStack != nullptr) && (originStack->GetNTrack() != 0u))
    {
      aStackedTrack = originStack->PopFromStack();
      if(targetStack != nullptr) { targetStack->PushToStack(aStackedTrack); }
      else            { urgentStack->PushToStack(aStackedTrack); }
    }
    else if(urgentStack->GetNTrack() != 0u)
    {
      aStackedTrack = urgentStack->PopFromStack();
      if(targetStack != nullptr) { targetStack->PushToStack(aStackedTrack); }
      else            { urgentStack->PushToStack(aStackedTrack); }
    }
  }
  return;
}

void G4StackManager::clear()
{
  ClearUrgentStack();
  ClearWaitingStack();
  for(G4int i=1; i<=numberOfAdditionalWaitingStacks; ++i)
  {
    ClearWaitingStack(i);
  }
}

void G4StackManager::ClearUrgentStack()
{
  urgentStack->clearAndDestroy();
}

void G4StackManager::ClearWaitingStack(G4int i)
{
  if(i==0)
  {
    waitingStack->clearAndDestroy();
  }
  else
  {
    if(i<=numberOfAdditionalWaitingStacks)
    {
      additionalWaitingStacks[i-1]->clearAndDestroy();
    }
  }
}

void G4StackManager::ClearPostponeStack()
{
  postponeStack->clearAndDestroy();
}

G4int G4StackManager::GetNTotalTrack() const
{
  std::size_t n = urgentStack->GetNTrack()
                + waitingStack->GetNTrack()
                + postponeStack->GetNTrack();
  for(G4int i=1; i<=numberOfAdditionalWaitingStacks; ++i)
  {
    n += additionalWaitingStacks[i-1]->GetNTrack();
  }
  return G4int(n);
}

G4int G4StackManager::GetNUrgentTrack() const
{
  return (G4int)urgentStack->GetNTrack();
}

G4int G4StackManager::GetNWaitingTrack(int i) const
{
  if(i==0)
  {
    return (G4int)waitingStack->GetNTrack();
  }
  
  if(i<=numberOfAdditionalWaitingStacks)
  {
    return (G4int)additionalWaitingStacks[i-1]->GetNTrack();
  }
 
  return 0;
}

G4int G4StackManager::GetNPostponedTrack() const
{
  return (G4int)postponeStack->GetNTrack();
}

void G4StackManager::SetVerboseLevel( G4int const value )
{
  verboseLevel = value;
}

void G4StackManager::SetUserStackingAction(G4UserStackingAction* value)
{
  userStackingAction = value;
  if(userStackingAction != nullptr)
  {
    userStackingAction->SetStackManager(this);
  }
}

G4ClassificationOfNewTrack G4StackManager::
DefaultClassification(G4Track* aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;
  if( aTrack->GetTrackStatus() == fPostponeToNextEvent )
  {
    classification = fPostpone;
  }
  return classification;
}
