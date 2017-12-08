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
//
// $Id: G4StackManager.cc 106992 2017-10-31 10:14:18Z gcosmo $
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
#ifdef G4_USESMARTSTACK
  urgentStack = new G4SmartTrackStack;
 // G4cout<<"+++ G4StackManager uses G4SmartTrackStack. +++"<<G4endl;
#else
  urgentStack = new G4TrackStack(5000);
//  G4cout<<"+++ G4StackManager uses ordinary G4TrackStack. +++"<<G4endl;
#endif
  waitingStack = new G4TrackStack(1000);
  postponeStack = new G4TrackStack(1000);
}

G4StackManager::~G4StackManager()
{
  if(userStackingAction) delete userStackingAction;

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

#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"

//Needed for temporal service
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

G4int G4StackManager::PushOneTrack(G4Track *newTrack,G4VTrajectory *newTrajectory)
{
  const G4ParticleDefinition* pd = newTrack->GetParticleDefinition();
  if(pd->GetParticleDefinitionID() < 0)
  {
#ifdef G4VERBOSE
    G4ExceptionDescription ED;
    if(verboseLevel>0) {
      ED << "A track without proper process manager is pushed into the track stack.\n"
         << " Particle name : " << pd->GetParticleName() << " -- ";
      if(newTrack->GetParentID()<0)
      { ED << "created by a primary particle generator."; }
      else
      { 
        const G4VProcess* vp = newTrack->GetCreatorProcess();
        if(vp)
        { ED << "created by " << vp->GetProcessName() << "."; }
        else
        { ED << "creaded by unknown process."; }
      }
    }
#endif
////  Temporal care of setting process manager for general ion.
    if(pd->IsGeneralIon())
    {
#ifdef G4VERBOSE
      if( verboseLevel > 0 ) {
        ED << "\n Process manager is temporally set, but this operation is thread-unsafe\n"
           << "and will be replaced with other methods at version 10.0.";
        G4Exception("G4StackManager::PushOneTrack","Event10051",JustWarning,ED);
      }
#endif
      G4ParticleDefinition* genericIon = G4ParticleTable::GetParticleTable()->GetGenericIon();
      G4ProcessManager* pman=0;
      if (genericIon!=0) pman = genericIon->GetProcessManager();
      if ((genericIon ==0) || (pman==0)){
        G4Exception( "G4IonTable::AddProcessManager()","PART10052", FatalException,
                   "Can not define process manager. GenericIon is not available.");
      }
      G4ParticleDefinition* ion = const_cast<G4ParticleDefinition*>(pd);
      ion->SetParticleDefinitionID(genericIon->GetParticleDefinitionID());
#ifdef G4VERBOSE
      if( verboseLevel > 1 )
      {
        G4ProcessManager* ionPman = ion->GetProcessManager();
        G4cout << "Now " << ion->GetParticleName() << " has a process manaegr at " << ionPman
               << " that is equivalent to " << pman << G4endl;
        G4ProcessVector* ionPvec = ionPman->GetProcessList();
        for(G4int ip1=0;ip1<ionPvec->size();ip1++)
        {
          G4cout << " " << ip1 << " - " << (*ionPvec)[ip1]->GetProcessName()
                 << " AtRest " << ionPman->GetAtRestIndex((*ionPvec)[ip1])
                 << ", AlongStep " << ionPman->GetAlongStepIndex((*ionPvec)[ip1])
                 << ", PostStep " << ionPman->GetPostStepIndex((*ionPvec)[ip1])
                 << G4endl;
        }
      }
#endif
    }
////  End of temporal care of setting process manager
#ifdef G4MUATOMS_INUSE
////  Setting process manager for muonic atom (a temporary solution)
    if(pd->IsMuonicAtom())
    {
#ifdef G4VERBOSE
      if( verboseLevel > 0 ) {
        ED << "\n Process manager is temporally set, but this operation is thread-unsafe\n"
           << "and will be replaced with other methods at version 10.0.";
        G4Exception("G4StackManager::PushOneTrack","Event10051",JustWarning,ED);
      }
#endif
      G4ParticleDefinition* genericMA = 
        G4ParticleTable::GetParticleTable()->GetGenericMuonicAtom();
      G4ProcessManager* pman=nullptr;
      if (genericMA!=nullptr) pman = genericMA->GetProcessManager();
      if ((genericMA == nullptr) || (pman== nullptr)){
        G4Exception( "G4IonTable::AddProcessManager()","PART10052", FatalException,
                   "Can not define process manager. GenericIon is not available.");
      }
      G4ParticleDefinition* muAtom = const_cast<G4ParticleDefinition*>(pd);
      muAtom->SetParticleDefinitionID(genericMA->GetParticleDefinitionID());
#ifdef G4VERBOSE
      if( verboseLevel > 1 )
      {
        G4ProcessManager* muAtomPman = muAtom->GetProcessManager();
        G4cout << "Now " << muAtom->GetParticleName() << " has a process manager at " << muAtomPman
               << " that is equivalent to " << pman << G4endl;
        G4ProcessVector* muAtomPvec = muAtomPman->GetProcessList();
        for(G4int ip1=0;ip1<muAtomPvec->size();ip1++)
        {
          G4cout << " " << ip1 << " - " << (*muAtomPvec)[ip1]->GetProcessName()
                 << " AtRest " << muAtomPman->GetAtRestIndex((*muAtomPvec)[ip1])
                 << ", AlongStep " << muAtomPman->GetAlongStepIndex((*muAtomPvec)[ip1])
                 << ", PostStep " << muAtomPman->GetPostStepIndex((*muAtomPvec)[ip1])
                 << G4endl;
        }
      }
#endif
    }
////  End of Setting process manager for muonic atom (a temporary solution)
#endif
    else
    {
#ifdef G4VERBOSE
      if( verboseLevel > 0 ) {
        ED << "\nThis track is deleted.";
        G4Exception("G4StackManager::PushOneTrack","Event10051",
                 JustWarning,ED);
      }
#endif
      delete newTrack;
      return GetNUrgentTrack();
    }
  }
    
  G4ClassificationOfNewTrack classification = DefaultClassification( newTrack ); 
  if(userStackingAction) 
  { classification = userStackingAction->ClassifyNewTrack( newTrack ); }

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
        if(i<1||i>numberOfAdditionalWaitingStacks) {
          G4ExceptionDescription ED;
          ED << "invalid classification " << classification << G4endl;
          G4Exception("G4StackManager::PushOneTrack","Event0051",
          FatalException,ED);
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
  if( verboseLevel > 1 )
  {
    G4cout << "### pop requested out of " 
         << GetNUrgentTrack() << " stacked tracks." << G4endl;
  }
#endif

  while( GetNUrgentTrack() == 0 )
  {
#ifdef G4VERBOSE
    if( verboseLevel > 1 ) G4cout << "### " << GetNWaitingTrack()
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
    if( verboseLevel > 1 ) G4cout << "     " << GetNUrgentTrack()
                      << " urgent tracks and " << GetNWaitingTrack()
                      << " waiting tracks." << G4endl;
#endif
    if( ( GetNUrgentTrack()==0 ) && ( GetNWaitingTrack()==0 ) ) return 0;
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
  
  if( !userStackingAction ) return;
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
        if(i<1||i>numberOfAdditionalWaitingStacks) {
          G4ExceptionDescription ED;
          ED << "invalid classification " << classification << G4endl;
          G4Exception("G4StackManager::ReClassify","Event0052",
                      FatalException,ED);
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
  
  urgentStack->clearAndDestroy(); // Set the urgentStack in a defined state. Not doing it would affect reproducibility.
  
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
      if(userStackingAction)
      { classification = userStackingAction->ClassifyNewTrack( aTrack ); }
      else
      { classification = DefaultClassification( aTrack ); }
      
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
            if(i<1||i>numberOfAdditionalWaitingStacks) {
              G4ExceptionDescription ED;
              ED << "invalid classification " << classification << G4endl;
              G4Exception("G4StackManager::PrepareNewEvent","Event0053",
                          FatalException,ED);
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
  if(origin==destination) return;
  if(origin==fKill) return;
  G4TrackStack* originStack = 0;
  switch(origin)
  {
    case fUrgent:
      originStack = 0;
      break;
    case fWaiting:
      originStack = waitingStack;
      break;
    case fPostpone:
      originStack = postponeStack;
      break;
    default:
      int i = origin - 10;
      if(i<=numberOfAdditionalWaitingStacks) originStack = additionalWaitingStacks[i-1];
      break;
  }
  
  if(destination==fKill)
  {
    if(originStack)
    { originStack->clearAndDestroy(); }
    else
    { urgentStack->clearAndDestroy(); }
  }
  else
  {
    G4TrackStack* targetStack = 0;
    switch(destination)
    {
      case fUrgent:
        targetStack = 0;
        break;
      case fWaiting:
        targetStack = waitingStack;
        break;
      case fPostpone:
        targetStack = postponeStack;
        break;
      default:
        int i = destination - 10;
        if(i<=numberOfAdditionalWaitingStacks) targetStack = additionalWaitingStacks[i-1];
        break;
    }
    if(originStack)
    {
      if(targetStack)
      { originStack->TransferTo(targetStack); }
      else
      { originStack->TransferTo(urgentStack); }
    }
    else
    { urgentStack->TransferTo(targetStack); }
  }
  return;
}

void G4StackManager::TransferOneStackedTrack(G4ClassificationOfNewTrack origin, G4ClassificationOfNewTrack destination)
{
  if(origin==destination) return;
  if(origin==fKill) return;
  G4TrackStack* originStack = 0;
  switch(origin)
  {
    case fUrgent:
      originStack = 0;
      break;
    case fWaiting:
      originStack = waitingStack;
      break;
    case fPostpone:
      originStack = postponeStack;
      break;
    default:
      int i = origin - 10;
      if(i<=numberOfAdditionalWaitingStacks) originStack = additionalWaitingStacks[i-1];
      break;
  }
  
  G4StackedTrack aStackedTrack;
  if(destination==fKill)
  {
    if( originStack && originStack->GetNTrack() ) {
      aStackedTrack = originStack->PopFromStack();
      delete aStackedTrack.GetTrack();
      delete aStackedTrack.GetTrajectory();
    }
    else if (urgentStack->GetNTrack() ) {
      aStackedTrack = urgentStack->PopFromStack();
      delete aStackedTrack.GetTrack();
      delete aStackedTrack.GetTrajectory();
    }
  }
  else
  {
    G4TrackStack* targetStack = 0;
    switch(destination)
    {
      case fUrgent:
        targetStack = 0;
        break;
      case fWaiting:
        targetStack = waitingStack;
        break;
      case fPostpone:
        targetStack = postponeStack;
        break;
      default:
        int i = destination - 10;
        if(i<=numberOfAdditionalWaitingStacks) targetStack = additionalWaitingStacks[i-1];
        break;
    }
    if(originStack && originStack->GetNTrack()) {
      aStackedTrack = originStack->PopFromStack();
      if(targetStack) { targetStack->PushToStack(aStackedTrack); }
      else            { urgentStack->PushToStack(aStackedTrack); }
    }
    else if(urgentStack->GetNTrack()) {
      aStackedTrack = urgentStack->PopFromStack();
      if(targetStack) { targetStack->PushToStack(aStackedTrack); }
      else            { urgentStack->PushToStack(aStackedTrack); }
    }
  }
  return;
}

void G4StackManager::clear()
{
  ClearUrgentStack();
  ClearWaitingStack();
  for(int i=1;i<=numberOfAdditionalWaitingStacks;i++) {ClearWaitingStack(i);}
}

void G4StackManager::ClearUrgentStack()
{
  urgentStack->clearAndDestroy();
}

void G4StackManager::ClearWaitingStack(int i)
{
  if(i==0) {
    waitingStack->clearAndDestroy();
  } else {
    if(i<=numberOfAdditionalWaitingStacks) additionalWaitingStacks[i-1]->clearAndDestroy();
  }
}

void G4StackManager::ClearPostponeStack()
{
  postponeStack->clearAndDestroy();
}

G4int G4StackManager::GetNTotalTrack() const
{
  int n = urgentStack->GetNTrack() + waitingStack->GetNTrack() + postponeStack->GetNTrack();
  for(int i=1;i<=numberOfAdditionalWaitingStacks;i++) {n += additionalWaitingStacks[i-1]->GetNTrack();}
  return n;
}

G4int G4StackManager::GetNUrgentTrack() const
{
  return urgentStack->GetNTrack();
}

G4int G4StackManager::GetNWaitingTrack(int i) const
{
  if(i==0) { return waitingStack->GetNTrack(); }
  else {
    if(i<=numberOfAdditionalWaitingStacks) { return additionalWaitingStacks[i-1]->GetNTrack();}
  }
  return 0;
}

G4int G4StackManager::GetNPostponedTrack() const
{
  return postponeStack->GetNTrack();
}

void G4StackManager::SetVerboseLevel( G4int const value )
{
  verboseLevel = value;
}

void G4StackManager::SetUserStackingAction(G4UserStackingAction* value)
{
	userStackingAction = value;
  if(userStackingAction) userStackingAction->SetStackManager(this);
}

G4ClassificationOfNewTrack G4StackManager::DefaultClassification(G4Track *aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;
  if( aTrack->GetTrackStatus() == fPostponeToNextEvent )
  { classification = fPostpone; }
  return classification;
}




