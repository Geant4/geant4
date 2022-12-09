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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITStepProcessor.hh"
#include "G4UImanager.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4ITTransportationManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4ParticleTable.hh"
#include "G4ITTrackingManager.hh"
#include "G4TrackingInformation.hh"
#include "G4IT.hh"
#include "G4ITNavigator.hh"             // Include from 'geometry'

#include "G4VITProcess.hh"
#include "G4VProcess.hh"
#include "G4ITTransportation.hh"

#include "G4ITTrackHolder.hh"

#include "G4ITSteppingVerbose.hh"
#include "G4VITSteppingVerbose.hh"
#include "G4ITTrackingInteractivity.hh"

#include "G4TrackingInformation.hh"

#include <iomanip>              // Include from 'system'
#include <vector>               // Include from 'system'

using namespace std;

static const std::size_t SizeOfSelectedDoItVector = 100;

template<typename T>
  inline bool IsInf(T value)
  {
    return std::numeric_limits<T>::has_infinity
        && value == std::numeric_limits<T>::infinity();
  }

//______________________________________________________________________________

G4ITStepProcessor::G4ITStepProcessor()
{
  fpVerbose = 0;
  // fpUserSteppingAction = 0 ;
  fStoreTrajectory = 0;
  fpTrackingManager = 0;
  fpNavigator = 0;
  kCarTolerance = -1.;
  fInitialized = false;
  fPreviousTimeStep = DBL_MAX;
  fILTimeStep = DBL_MAX;
  fpTrackContainer = 0;

  CleanProcessor();
  ResetSecondaries();
}

//______________________________________________________________________________

//G4ITStepProcessor::
G4ITStepProcessorState::G4ITStepProcessorState() :
    G4ITStepProcessorState_Lock(),
    fSelectedAtRestDoItVector(G4VITProcess::GetMaxProcessIndex(), 0),
    fSelectedPostStepDoItVector(G4VITProcess::GetMaxProcessIndex(), 0)
{
  fPhysicalStep = -1.;
  fPreviousStepSize = -1.;

  fSafety = -1.;
  fProposedSafety = -1.;
  fEndpointSafety = -1;

  fStepStatus = fUndefined;

  fTouchableHandle = 0;
}

//______________________________________________________________________________

//G4ITStepProcessor::
G4ITStepProcessorState::G4ITStepProcessorState(const G4ITStepProcessorState& right) :
    G4ITStepProcessorState_Lock(),
    fSelectedAtRestDoItVector(right.fSelectedAtRestDoItVector),
    fSelectedPostStepDoItVector(right.fSelectedPostStepDoItVector)
{
  fPhysicalStep = right.fPhysicalStep;
  fPreviousStepSize = right.fPreviousStepSize;

  fSafety = right.fSafety;
  fProposedSafety = right.fProposedSafety;
  fEndpointSafety = right.fEndpointSafety;

  fStepStatus = right.fStepStatus;

  fTouchableHandle = right.fTouchableHandle;
}

//______________________________________________________________________________

//G4ITStepProcessor::
G4ITStepProcessorState&
//G4ITStepProcessor::
G4ITStepProcessorState::operator=(const G4ITStepProcessorState& right)
{
  if(this == &right) return *this;

  fSelectedAtRestDoItVector.clear();
  fSelectedAtRestDoItVector = right.fSelectedAtRestDoItVector;
  fSelectedPostStepDoItVector.clear();
  fSelectedPostStepDoItVector = right.fSelectedPostStepDoItVector;

  fPhysicalStep = right.fPhysicalStep;
  fPreviousStepSize = right.fPreviousStepSize;

  fSafety = right.fSafety;
  fProposedSafety = right.fProposedSafety;
  fEndpointSafety = right.fEndpointSafety;

  fStepStatus = right.fStepStatus;

  fTouchableHandle = right.fTouchableHandle;
  return *this;
}

//______________________________________________________________________________

//G4ITStepProcessor::
G4ITStepProcessorState::~G4ITStepProcessorState()
{
  ;
}

//______________________________________________________________________________

void G4ITStepProcessor::ClearProcessInfo()
{
  std::map<const G4ParticleDefinition*, ProcessGeneralInfo*>::iterator it;

  for(it = fProcessGeneralInfoMap.begin(); it != fProcessGeneralInfoMap.end();
      it++)
  {
    if(it->second)
    {
      delete it->second;
      it->second = 0;
    }
  }

  fProcessGeneralInfoMap.clear();
}

//______________________________________________________________________________

void G4ITStepProcessor::ForceReInitialization()
{
  fInitialized = false;
  ClearProcessInfo();
  Initialize();
}

//______________________________________________________________________________

void G4ITStepProcessor::Initialize()
{
  CleanProcessor();
  if(fInitialized) return;
  //    ActiveOnlyITProcess();

  SetNavigator(G4ITTransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking());

  fPhysIntLength = DBL_MAX;
  kCarTolerance = 0.5
      * G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if(fpVerbose == 0)
  {
    G4ITTrackingInteractivity* interactivity = fpTrackingManager->GetInteractivity();

    if(interactivity)
    {
      fpVerbose = interactivity->GetSteppingVerbose();
      fpVerbose->SetStepProcessor(this);
    }
  }

  fpTrackContainer = G4ITTrackHolder::Instance();

  fInitialized = true;
}
//______________________________________________________________________________

G4ITStepProcessor::~G4ITStepProcessor()
{
  if(fpStep)
  {
    fpStep->DeleteSecondaryVector();
    delete fpStep;
  }

  if(fpSecondary) delete fpSecondary;
  ClearProcessInfo();
  //G4ITTransportationManager::DeleteInstance();

  //    if(fpUserSteppingAction)             delete fpUserSteppingAction;
}
//______________________________________________________________________________
// should not be used
G4ITStepProcessor::G4ITStepProcessor(const G4ITStepProcessor& rhs)
{
  fpVerbose = rhs.fpVerbose;
  fStoreTrajectory = rhs.fStoreTrajectory;

  //    fpUserSteppingAction = 0 ;
  fpTrackingManager = 0;
  fpNavigator = 0;
  fInitialized = false;

  kCarTolerance = rhs.kCarTolerance;
  fInitialized = false;
  fPreviousTimeStep = DBL_MAX;

  CleanProcessor();
  ResetSecondaries();
  fpTrackContainer = 0;
  fILTimeStep = DBL_MAX;
}

//______________________________________________________________________________

void G4ITStepProcessor::ResetLeadingTracks()
{
  fLeadingTracks.Reset();
}

//______________________________________________________________________________

void G4ITStepProcessor::PrepareLeadingTracks()
{
  fLeadingTracks.PrepareLeadingTracks();
}

//______________________________________________________________________________

G4ITStepProcessor& G4ITStepProcessor::operator=(const G4ITStepProcessor& rhs)
{
  if(this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

//______________________________________________________________________________

void G4ITStepProcessor::ActiveOnlyITProcess()
{
  // Method not used for the time being
#ifdef debug
  G4cout<<"G4ITStepProcessor::CloneProcesses: is called"<<G4endl;
#endif

  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable
      ->GetIterator();

  theParticleIterator->reset();
  // TODO : Ne faire la boucle que sur les IT **** !!!
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pm = particle->GetProcessManager();

    if(!pm)
    {
      G4cerr << "ERROR - G4ITStepProcessor::GetProcessNumber()" << G4endl<< "        ProcessManager is NULL for particle = "
      << particle->GetParticleName() << ", PDG_code = "
      << particle->GetPDGEncoding() << G4endl;
      G4Exception("G4ITStepProcessor::GetProcessNumber()", "ITStepProcessor0001",
          FatalException, "Process Manager is not found.");
      return;
    }

    ActiveOnlyITProcess(pm);
  }
}

//______________________________________________________________________________

void G4ITStepProcessor::ActiveOnlyITProcess(G4ProcessManager* processManager)
{
  // Method not used for the time being
  G4ProcessVector* processVector = processManager->GetProcessList();

  G4VITProcess* itProcess = 0;
  for(G4int i = 0; i < (G4int)processVector->size(); ++i)
  {
    G4VProcess* base_process = (*processVector)[i];
    itProcess = dynamic_cast<G4VITProcess*>(base_process);

    if(!itProcess)
    {
      processManager->SetProcessActivation(base_process, false);
    }
  }
}

//______________________________________________________________________________

void G4ITStepProcessor::SetupGeneralProcessInfo(G4ParticleDefinition* particle,
                                                G4ProcessManager* pm)
{

#ifdef debug
  G4cout<<"G4ITStepProcessor::GetProcessNumber: is called track"<<G4endl;
#endif
  if(!pm)
  {
    G4cerr << "ERROR - G4SteppingManager::GetProcessNumber()" << G4endl<< "        ProcessManager is NULL for particle = "
    << particle->GetParticleName() << ", PDG_code = "
    << particle->GetPDGEncoding() << G4endl;
    G4Exception("G4SteppingManager::GetProcessNumber()", "ITStepProcessor0002",
        FatalException, "Process Manager is not found.");
    return;
  }

  std::map<const G4ParticleDefinition*, ProcessGeneralInfo*>::iterator it =
  fProcessGeneralInfoMap.find(particle);
  if(it != fProcessGeneralInfoMap.end())
  {
    G4Exception("G4SteppingManager::SetupGeneralProcessInfo()",
        "ITStepProcessor0003",
        FatalException, "Process info already registered.");
    return;
  }

  // here used as temporary
  fpProcessInfo = new ProcessGeneralInfo();

  // AtRestDoits
  fpProcessInfo->MAXofAtRestLoops = pm->GetAtRestProcessVector()->entries();
  fpProcessInfo->fpAtRestDoItVector = pm->GetAtRestProcessVector(typeDoIt);
  fpProcessInfo->fpAtRestGetPhysIntVector =
  pm->GetAtRestProcessVector(typeGPIL);
#ifdef debug
  G4cout << "G4ITStepProcessor::GetProcessNumber: #ofAtRest="
  << fpProcessInfo->MAXofAtRestLoops << G4endl;
#endif

  // AlongStepDoits
  fpProcessInfo->MAXofAlongStepLoops =
  pm->GetAlongStepProcessVector()->entries();
  fpProcessInfo->fpAlongStepDoItVector =
  pm->GetAlongStepProcessVector(typeDoIt);
  fpProcessInfo->fpAlongStepGetPhysIntVector =
  pm->GetAlongStepProcessVector(typeGPIL);
#ifdef debug
  G4cout << "G4ITStepProcessor::GetProcessNumber:#ofAlongStp="
  << fpProcessInfo->MAXofAlongStepLoops << G4endl;
#endif

  // PostStepDoits
  fpProcessInfo->MAXofPostStepLoops =
  pm->GetPostStepProcessVector()->entries();
  fpProcessInfo->fpPostStepDoItVector = pm->GetPostStepProcessVector(typeDoIt);
  fpProcessInfo->fpPostStepGetPhysIntVector =
  pm->GetPostStepProcessVector(typeGPIL);
#ifdef debug
  G4cout << "G4ITStepProcessor::GetProcessNumber: #ofPostStep="
  << fpProcessInfo->MAXofPostStepLoops << G4endl;
#endif

  if (SizeOfSelectedDoItVector<fpProcessInfo->MAXofAtRestLoops ||
      SizeOfSelectedDoItVector<fpProcessInfo->MAXofAlongStepLoops ||
      SizeOfSelectedDoItVector<fpProcessInfo->MAXofPostStepLoops )
  {
    G4cerr << "ERROR - G4ITStepProcessor::GetProcessNumber()" << G4endl
    << "        SizeOfSelectedDoItVector= " << SizeOfSelectedDoItVector
    << " ; is smaller then one of MAXofAtRestLoops= "
    << fpProcessInfo->MAXofAtRestLoops << G4endl
    << "        or MAXofAlongStepLoops= " << fpProcessInfo->MAXofAlongStepLoops
    << " or MAXofPostStepLoops= " << fpProcessInfo->MAXofPostStepLoops << G4endl;
    G4Exception("G4ITStepProcessor::GetProcessNumber()",
        "ITStepProcessor0004", FatalException,
        "The array size is smaller than the actual No of processes.");
  }

  if(!fpProcessInfo->fpAtRestDoItVector &&
      !fpProcessInfo->fpAlongStepDoItVector &&
      !fpProcessInfo->fpPostStepDoItVector)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No DoIt process found ";
    G4Exception("G4ITStepProcessor::DoStepping","ITStepProcessor0005",
        FatalErrorInArgument,exceptionDescription);
    return;
  }

  if(fpProcessInfo->fpAlongStepGetPhysIntVector
      && fpProcessInfo->MAXofAlongStepLoops>0)
  {
    fpProcessInfo->fpTransportation = dynamic_cast<G4ITTransportation*>
    ((*fpProcessInfo->fpAlongStepGetPhysIntVector)
        [G4int(fpProcessInfo->MAXofAlongStepLoops-1)]);

    if(fpProcessInfo->fpTransportation == 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No transportation process found ";
      G4Exception("G4ITStepProcessor::SetupGeneralProcessInfo",
          "ITStepProcessor0006",
          FatalErrorInArgument,exceptionDescription);
    }
  }
  fProcessGeneralInfoMap[particle] = fpProcessInfo;
  //    fpProcessInfo = 0;
}

//______________________________________________________________________________

void G4ITStepProcessor::SetTrack(G4Track* track)
{
  fpTrack = track;
  if(fpTrack)
  {
    fpITrack = GetIT(fpTrack);
    fpStep = const_cast<G4Step*>(fpTrack->GetStep());

    if(fpITrack)
    {
      fpTrackingInfo = fpITrack->GetTrackingInfo();
    }
    else
    {
      fpTrackingInfo = 0;
      G4cerr << "Track ID : " << fpTrack->GetTrackID() << G4endl;

      G4ExceptionDescription errMsg;
      errMsg << "No IT pointer was attached to the track you try to process.";
      G4Exception("G4ITStepProcessor::SetTrack",
                  "ITStepProcessor0007",
                  FatalErrorInArgument,
                  errMsg);
    }
  }
  else
  {
    fpITrack = 0;
    fpStep = 0;
  }
}
//______________________________________________________________________________

void G4ITStepProcessor::GetProcessInfo()
{
  G4ParticleDefinition* particle = fpTrack->GetDefinition();
  std::map<const G4ParticleDefinition*, ProcessGeneralInfo*>::iterator it =
      fProcessGeneralInfoMap.find(particle);

  if(it == fProcessGeneralInfoMap.end())
  {
    SetupGeneralProcessInfo(particle,
                            fpTrack->GetDefinition()->GetProcessManager());
    if(fpProcessInfo == 0)
    {
      G4ExceptionDescription exceptionDescription("...");
      G4Exception("G4ITStepProcessor::GetProcessNumber",
                  "ITStepProcessor0008",
                  FatalErrorInArgument,
                  exceptionDescription);
      return;
    }
  }
  else
  {
    fpProcessInfo = it->second;
  }
}

//______________________________________________________________________________

void G4ITStepProcessor::SetupMembers()
{
  fpSecondary = fpStep->GetfSecondary();
  fpPreStepPoint = fpStep->GetPreStepPoint();
  fpPostStepPoint = fpStep->GetPostStepPoint();

  fpState = (G4ITStepProcessorState*) fpITrack->GetTrackingInfo()
      ->GetStepProcessorState();

  GetProcessInfo();
  ResetSecondaries();
}

//______________________________________________________________________________

void G4ITStepProcessor::ResetSecondaries()
{
  // Reset the secondary particles
  fN2ndariesAtRestDoIt = 0;
  fN2ndariesAlongStepDoIt = 0;
  fN2ndariesPostStepDoIt = 0;
}

//______________________________________________________________________________

void G4ITStepProcessor::GetAtRestIL()
{
  // Select the rest process which has the shortest time before
  // it is invoked. In rest processes, GPIL()
  // returns the time before a process occurs.
  G4double lifeTime(DBL_MAX), shortestLifeTime (DBL_MAX);

  fAtRestDoItProcTriggered = 0;
  shortestLifeTime = DBL_MAX;

  unsigned int NofInactiveProc=0;

  for( G4int ri=0; ri < (G4int)fpProcessInfo->MAXofAtRestLoops; ++ri )
  {
    fpCurrentProcess = dynamic_cast<G4VITProcess*>((*fpProcessInfo->fpAtRestGetPhysIntVector)[ri]);
    if (fpCurrentProcess== 0)
    {
      (fpState->fSelectedAtRestDoItVector)[ri] = InActivated;
      NofInactiveProc++;
      continue;
    } // NULL means the process is inactivated by a user on fly.

    fCondition=NotForced;
    fpCurrentProcess->SetProcessState(
        fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));

    lifeTime = fpCurrentProcess->AtRestGPIL( *fpTrack, &fCondition );
    fpCurrentProcess->ResetProcessState();

    if(fCondition==Forced)
    {
      (fpState->fSelectedAtRestDoItVector)[ri] = Forced;
    }
    else
    {
      (fpState->fSelectedAtRestDoItVector)[ri] = InActivated;
      if(lifeTime < shortestLifeTime )
      {
        shortestLifeTime = lifeTime;
        fAtRestDoItProcTriggered = G4int(ri);
      }
    }
  }

  (fpState->fSelectedAtRestDoItVector)[fAtRestDoItProcTriggered] = NotForced;

//  G4cout << " --> Selected at rest process : "
//         << (*fpProcessInfo->fpAtRestGetPhysIntVector)[fAtRestDoItProcTriggered]
//            ->GetProcessName()
//         << G4endl;

  fTimeStep = shortestLifeTime;

  // at least one process is necessary to destroy the particle
  // exit with warning
  if(NofInactiveProc==fpProcessInfo->MAXofAtRestLoops)
  {
    G4cerr << "ERROR - G4ITStepProcessor::InvokeAtRestDoItProcs()" << G4endl
    << "        No AtRestDoIt process is active!" << G4endl;
  }
}

//_________________________________________________________________________
G4double G4ITStepProcessor::ComputeInteractionLength(double previousTimeStep)
{
  G4TrackManyList* mainList = fpTrackContainer->GetMainList();
  G4TrackManyList::iterator it = mainList ->begin();
  G4TrackManyList::iterator end = mainList ->end();

  SetPreviousStepTime(previousTimeStep);

  fILTimeStep = DBL_MAX;

  for (; it != end; )
  {
    G4Track * track = *it;

#ifdef DEBUG
    G4cout << "*CIL* " << GetIT(track)->GetName()
        << " ID: " << track->GetTrackID()
        <<  " at time : " << track->GetGlobalTime()
        << G4endl;
#endif

    ++it;
    DefinePhysicalStepLength(track);

    ExtractILData();
  }

  return fILTimeStep;
}

//_________________________________________________________________________

void G4ITStepProcessor::ExtractILData()
{
  assert(fpTrack != 0);
  if (fpTrack == 0)
  {
    CleanProcessor();
    return;
  }

  // assert(fpTrack->GetTrackStatus() != fStopAndKill);

  if (fpTrack->GetTrackStatus() == fStopAndKill)
  {
//    trackContainer->GetMainList()->pop(fpTrack);
    fpTrackingManager->EndTracking(fpTrack);
    CleanProcessor();
    return;
  }

  if (IsInf(fTimeStep))
  {
    // G4cout << "!!!!!!!!!!!! IS INF " << track->GetTrackID() << G4endl;
    CleanProcessor();
    return;
  }
  else if (fTimeStep < fILTimeStep - DBL_EPSILON)
  {
    // G4cout << "!!!!!!!!!!!! TEMPS DIFFERENTS "
    //    << track->GetTrackID() << G4endl;

    fLeadingTracks.Reset();

    fILTimeStep = GetInteractionTime();

//    G4cout << "Will set leading step to true for time  :"
//           << SP -> GetInteractionTime() << " against fTimeStep : "
//           << G4BestUnit(fILTimeStep, "Time") << " the trackID is : " << track->GetTrackID()
//           << G4endl;

//    GetIT(fpTrack)->GetTrackingInfo()->SetLeadingStep(true);
    fLeadingTracks.Push(fpTrack);
  }
  else if(fabs(fILTimeStep - fTimeStep) < DBL_EPSILON )
  {

    // G4cout << "!!!!!!!!!!!! MEME TEMPS " << track->GetTrackID() << G4endl;
    // G4cout << "Will set leading step to true for time  :"
    //        << SP -> GetInteractionTime() << " against fTimeStep : "
    //        << fTimeStep << " the trackID is : " << track->GetTrackID()<< G4endl;
//    GetIT(fpTrack)->GetTrackingInfo()->SetLeadingStep(true);
    fLeadingTracks.Push(fpTrack);
  }
  // else
  // {
  //     G4cout << "!!!! Bigger time : " << "currentTime : "<<fILTimeStep
  //  << " proposedTime : " << SP -> GetInteractionTime() << G4endl;
  // }

  CleanProcessor();
}

//___________________________________________________________________________

void G4ITStepProcessor::DefinePhysicalStepLength(G4Track* track)
{
  SetTrack(track);
  DoDefinePhysicalStepLength();
}

//______________________________________________________________________________

void G4ITStepProcessor::SetInitialStep()
{
  // DEBUG
  //    G4cout << "SetInitialStep for : " << fpITrack-> GetName() << G4endl;
  //________________________________________________________
  // Initialize geometry

  if(!fpTrack->GetTouchableHandle())
  {
    //==========================================================================
    // Create navigator state and Locate particle in geometry
    //==========================================================================
    /*
     fpNavigator->NewNavigatorStateAndLocate(fpTrack->GetPosition(),
     fpTrack->GetMomentumDirection());

     fpITrack->GetTrackingInfo()->
     SetNavigatorState(fpNavigator->GetNavigatorState());
     */
    fpNavigator->NewNavigatorState();
    fpITrack->GetTrackingInfo()->SetNavigatorState(fpNavigator
        ->GetNavigatorState());

    G4ThreeVector direction = fpTrack->GetMomentumDirection();
    fpNavigator->LocateGlobalPointAndSetup(fpTrack->GetPosition(),
                                           &direction,
                                           false,
                                           false); // was false, false

    fpState->fTouchableHandle = fpNavigator->CreateTouchableHistory();

    fpTrack->SetTouchableHandle(fpState->fTouchableHandle);
    fpTrack->SetNextTouchableHandle(fpState->fTouchableHandle);
  }
  else
  {
    fpState->fTouchableHandle = fpTrack->GetTouchableHandle();
    fpTrack->SetNextTouchableHandle(fpState->fTouchableHandle);

    //==========================================================================
    // Create OR set navigator state
    //==========================================================================

    if(fpITrack->GetTrackingInfo()->GetNavigatorState())
    {
      fpNavigator->SetNavigatorState(fpITrack->GetTrackingInfo()
          ->GetNavigatorState());
      fpITrack->GetTrackingInfo()->SetNavigatorState(fpNavigator
          ->GetNavigatorState());
    }
    else
    {
      fpNavigator->NewNavigatorState(*((G4TouchableHistory*) fpState
          ->fTouchableHandle()));
      fpITrack->GetTrackingInfo()->SetNavigatorState(fpNavigator
          ->GetNavigatorState());
    }

    G4VPhysicalVolume* oldTopVolume =
        fpTrack->GetTouchableHandle()->GetVolume();

    //==========================================================================
    // Locate particle in geometry
    //==========================================================================

//    G4VPhysicalVolume* newTopVolume =
//        fpNavigator->LocateGlobalPointAndSetup(
//            fpTrack->GetPosition(),
//            &fpTrack->GetMomentumDirection(),
//            true, false);

    G4VPhysicalVolume* newTopVolume =
        fpNavigator->ResetHierarchyAndLocate(fpTrack->GetPosition(),
                                             fpTrack->GetMomentumDirection(),
                                             *((G4TouchableHistory*) fpTrack
                                                 ->GetTouchableHandle()()));

    if(newTopVolume != oldTopVolume || oldTopVolume->GetRegularStructureId()
        == 1)
    {
      fpState->fTouchableHandle = fpNavigator->CreateTouchableHistory();
      fpTrack->SetTouchableHandle(fpState->fTouchableHandle);
      fpTrack->SetNextTouchableHandle(fpState->fTouchableHandle);
    }
  }

  fpCurrentVolume = fpState->fTouchableHandle->GetVolume();

  //________________________________________________________
  // If the primary track has 'Suspend' or 'PostponeToNextEvent' state,
  // set the track state to 'Alive'.
  if((fpTrack->GetTrackStatus() == fSuspend) || (fpTrack->GetTrackStatus()
      == fPostponeToNextEvent))
  {
    fpTrack->SetTrackStatus(fAlive);
  }

  //HoangTRAN: it's better to check the status here
  if(fpTrack->GetTrackStatus() == fStopAndKill) return;

  // If the primary track has 'zero' kinetic energy, set the track
  // state to 'StopButAlive'.
  if(fpTrack->GetKineticEnergy() <= 0.0)
  {
    fpTrack->SetTrackStatus(fStopButAlive);
  }
  //________________________________________________________
  // Set vertex information of G4Track at here
  if(fpTrack->GetCurrentStepNumber() == 0)
  {
    fpTrack->SetVertexPosition(fpTrack->GetPosition());
    fpTrack->SetVertexMomentumDirection(fpTrack->GetMomentumDirection());
    fpTrack->SetVertexKineticEnergy(fpTrack->GetKineticEnergy());
    fpTrack->SetLogicalVolumeAtVertex(fpTrack->GetVolume()->GetLogicalVolume());
  }
  //________________________________________________________
  // If track is already outside the world boundary, kill it
  if(fpCurrentVolume == 0)
  {
    // If the track is a primary, stop processing
    if(fpTrack->GetParentID() == 0)
    {
      G4cerr << "ERROR - G4ITStepProcessor::SetInitialStep()" << G4endl<< "        Primary particle starting at - "
      << fpTrack->GetPosition()
      << " - is outside of the world volume." << G4endl;
      G4Exception("G4ITStepProcessor::SetInitialStep()", "ITStepProcessor0011",
          FatalException, "Primary vertex outside of the world!");
    }

    fpTrack->SetTrackStatus( fStopAndKill );
    G4cout << "WARNING - G4ITStepProcessor::SetInitialStep()" << G4endl
    << "          Initial track position is outside world! - "
    << fpTrack->GetPosition() << G4endl;
  }
  else
  {
    // Initial set up for attribues of 'Step'
    fpStep->InitializeStep( fpTrack );
  }

  fpState->fStepStatus = fUndefined;
}
//______________________________________________________________________________

void G4ITStepProcessor::InitDefineStep()
{

  if(!fpStep)
  {
    // Create new Step and give it to the track
    fpStep = new G4Step();
    fpTrack->SetStep(fpStep);
    fpSecondary = fpStep->NewSecondaryVector();

    // Create new state and set it in the trackingInfo
    fpState = new G4ITStepProcessorState();
    fpITrack->GetTrackingInfo()->SetStepProcessorState((G4ITStepProcessorState_Lock*) fpState);

    SetupMembers();
    SetInitialStep();

    fpTrackingManager->StartTracking(fpTrack);
  }
  else
  {
    SetupMembers();

    fpState->fPreviousStepSize = fpTrack->GetStepLength();
    /***
     // Send G4Step information to Hit/Dig if the volume is sensitive
     fpCurrentVolume = fpStep->GetPreStepPoint()->GetPhysicalVolume();
     StepControlFlag =  fpStep->GetControlFlag();
     if( fpCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation)
     {
     fpSensitive = fpStep->GetPreStepPoint()->
     GetSensitiveDetector();

     //            if( fSensitive != 0 ) {
     //              fSensitive->Hit(fStep);
     //            }
     }
     ***/
    // Store last PostStepPoint to PreStepPoint, and swap current and next
    // volume information of G4Track. Reset total energy deposit in one Step.
    fpStep->CopyPostToPreStepPoint();
    fpStep->ResetTotalEnergyDeposit();

    //JA Set the volume before it is used (in DefineStepLength() for User Limit)
    fpCurrentVolume = fpStep->GetPreStepPoint()->GetPhysicalVolume();
    /*
     G4cout << G4endl;
     G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
     G4cout << "PreStepPoint Volume : "
     << fpCurrentVolume->GetName() << G4endl;
     G4cout << "Track Touchable : "
     << fpTrack->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
     G4cout << "Track NextTouchable : "
     << fpTrack->GetNextTouchableHandle()->GetVolume()->GetName()
     << G4endl;
     */
    // Reset the step's auxiliary points vector pointer
    fpStep->SetPointerToVectorOfAuxiliaryPoints(0);

    // Switch next touchable in track to current one
    fpTrack->SetTouchableHandle(fpTrack->GetNextTouchableHandle());
    fpState->fTouchableHandle = fpTrack->GetTouchableHandle();
    fpTrack->SetNextTouchableHandle(fpState->fTouchableHandle);

    //! ADDED BACK
    /*
     G4VPhysicalVolume* oldTopVolume =
     fpTrack->GetTouchableHandle()->GetVolume();
     fpNavigator->SetNavigatorState(
     fpITrack->GetTrackingInfo()->GetNavigatorState());

     G4VPhysicalVolume* newTopVolume = fpNavigator->ResetHierarchyAndLocate(
     fpTrack->GetPosition(), fpTrack->GetMomentumDirection(),
     *((G4TouchableHistory*) fpTrack->GetTouchableHandle()()));

     //  G4VPhysicalVolume* newTopVolume=
     //     fpNavigator->LocateGlobalPointAndSetup(fpTrack->GetPosition(),
     //                                            &fpTrack->GetMomentumDirection(),
     //                                            true, false);

     //        G4cout << "New Top Volume : " << newTopVolume->GetName() << G4endl;

     if (newTopVolume != oldTopVolume || oldTopVolume->GetRegularStructureId()
     == 1)
     {
     fpState->fTouchableHandle = fpNavigator->CreateTouchableHistory();
     fpTrack->SetTouchableHandle(fpState->fTouchableHandle);
     fpTrack->SetNextTouchableHandle(fpState->fTouchableHandle);
     }

     */
    //! ADDED BACK
    //==========================================================================
    // Only reset navigator state + reset volume hierarchy (internal)
    // No need to relocate
    //==========================================================================
    fpNavigator->SetNavigatorState(fpITrack->GetTrackingInfo()
        ->GetNavigatorState());
  }
}

//______________________________________________________________________________

// ------------------------------------------------------------------------
//	Compute Interaction Length
// ------------------------------------------------------------------------
void G4ITStepProcessor::DoDefinePhysicalStepLength()
{

  InitDefineStep();

#ifdef G4VERBOSE
  // !!!!! Verbose
  if(fpVerbose) fpVerbose->DPSLStarted();
#endif

  G4TrackStatus trackStatus = fpTrack->GetTrackStatus();

  if(trackStatus == fStopAndKill)
  {
    return;
  }

  if(trackStatus == fStopButAlive)
  {
    fpITrack->GetTrackingInfo()->SetNavigatorState(fpNavigator
        ->GetNavigatorState());
    fpNavigator->ResetNavigatorState();
    return GetAtRestIL();
  }

  // Find minimum Step length and corresponding time
  // demanded by active disc./cont. processes

  // ReSet the counter etc.
  fpState->fPhysicalStep = DBL_MAX; // Initialize by a huge number
  fPhysIntLength = DBL_MAX; // Initialize by a huge number

  G4double proposedTimeStep = DBL_MAX;
  G4VProcess* processWithPostStepGivenByTimeStep(0);

  // GPIL for PostStep
  fPostStepDoItProcTriggered = fpProcessInfo->MAXofPostStepLoops;
  fPostStepAtTimeDoItProcTriggered = fpProcessInfo->MAXofPostStepLoops;

  //  G4cout << "fpProcessInfo->MAXofPostStepLoops : "
  //         << fpProcessInfo->MAXofPostStepLoops
  //         << " mol : " << fpITrack -> GetName()
  //         << " id : " << fpTrack->GetTrackID()
  //         << G4endl;

  for(std::size_t np = 0; np < fpProcessInfo->MAXofPostStepLoops; ++np)
  {
    fpCurrentProcess = dynamic_cast<G4VITProcess*>((*fpProcessInfo
        ->fpPostStepGetPhysIntVector)[(G4int)np]);
    if(fpCurrentProcess == 0)
    {
      (fpState->fSelectedPostStepDoItVector)[np] = InActivated;
      continue;
    } // NULL means the process is inactivated by a user on fly.

    fCondition = NotForced;
    fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess
        ->GetProcessID()));

    // G4cout << "Is going to call : "
    //        << fpCurrentProcess -> GetProcessName()
    //        << G4endl;
    fPhysIntLength = fpCurrentProcess->PostStepGPIL(*fpTrack,
                                                    fpState->fPreviousStepSize,
                                                    &fCondition);

#ifdef G4VERBOSE
    // !!!!! Verbose
    if(fpVerbose) fpVerbose->DPSLPostStep();
#endif

    fpCurrentProcess->ResetProcessState();
    //fpCurrentProcess->SetProcessState(0);

    switch(fCondition)
    {
      case ExclusivelyForced: // Will need special treatment
        (fpState->fSelectedPostStepDoItVector)[np] = ExclusivelyForced;
        fpState->fStepStatus = fExclusivelyForcedProc;
        fpStep->GetPostStepPoint()->SetProcessDefinedStep(fpCurrentProcess);
        break;

      case Conditionally:
        //	     (fpState->fSelectedPostStepDoItVector)[np] = Conditionally;
        G4Exception("G4ITStepProcessor::DefinePhysicalStepLength()",
                    "ITStepProcessor0008",
                    FatalException,
                    "This feature is no more supported");
        break;

      case Forced:
        (fpState->fSelectedPostStepDoItVector)[np] = Forced;
        break;

      case StronglyForced:
        (fpState->fSelectedPostStepDoItVector)[np] = StronglyForced;
        break;

      default:
        (fpState->fSelectedPostStepDoItVector)[np] = InActivated;
        break;
    }

    if(fCondition == ExclusivelyForced)
    {
      for(std::size_t nrest = np + 1; nrest < fpProcessInfo->MAXofPostStepLoops;
          ++nrest)
      {
        (fpState->fSelectedPostStepDoItVector)[nrest] = InActivated;
      }
      return; // Please note the 'return' at here !!!
    }
    else
    {
      if(fPhysIntLength < fpState->fPhysicalStep)
      {
        // To avoid checking whether the process is actually
        // proposing a time step, the returned time steps are
        // negative (just for tagging)
        if(fpCurrentProcess->ProposesTimeStep())
        {
          fPhysIntLength *= -1;
          if(fPhysIntLength < proposedTimeStep)
          {
            proposedTimeStep = fPhysIntLength;
            fPostStepAtTimeDoItProcTriggered = np;
            processWithPostStepGivenByTimeStep = fpCurrentProcess;
          }
        }
        else
        {
          fpState->fPhysicalStep = fPhysIntLength;
          fpState->fStepStatus = fPostStepDoItProc;
          fPostStepDoItProcTriggered = G4int(np);
          fpStep->GetPostStepPoint()->SetProcessDefinedStep(fpCurrentProcess);
        }
      }
    }
  }

  // GPIL for AlongStep
  fpState->fProposedSafety = DBL_MAX;
  G4double safetyProposedToAndByProcess = fpState->fProposedSafety;

  for(std::size_t kp = 0; kp < fpProcessInfo->MAXofAlongStepLoops; ++kp)
  {
    fpCurrentProcess = dynamic_cast<G4VITProcess*>((*fpProcessInfo
        ->fpAlongStepGetPhysIntVector)[(G4int)kp]);
    if(fpCurrentProcess == 0) continue;
    // NULL means the process is inactivated by a user on fly.

    fpCurrentProcess->SetProcessState(
        fpTrackingInfo->GetProcessState(fpCurrentProcess->GetProcessID()));
    fPhysIntLength =
        fpCurrentProcess->AlongStepGPIL(*fpTrack,
                                        fpState->fPreviousStepSize,
                                        fpState->fPhysicalStep,
                                        safetyProposedToAndByProcess,
                                        &fGPILSelection);

#ifdef G4VERBOSE
    // !!!!! Verbose
    if(fpVerbose) fpVerbose->DPSLAlongStep();
#endif

    if(fPhysIntLength < fpState->fPhysicalStep)
    {
      fpState->fPhysicalStep = fPhysIntLength;
      // Should save PS and TS in IT

      // Check if the process wants to be the GPIL winner. For example,
      // multi-scattering proposes Step limit, but won't be the winner.
      if(fGPILSelection == CandidateForSelection)
      {
        fpState->fStepStatus = fAlongStepDoItProc;
        fpStep->GetPostStepPoint()->SetProcessDefinedStep(fpCurrentProcess);
      }

      // Transportation is assumed to be the last process in the vector
      if(kp == fpProcessInfo->MAXofAlongStepLoops - 1)
      {
        fpTransportation = dynamic_cast<G4ITTransportation*>(fpCurrentProcess);

        if(!fpTransportation)
        {
          G4ExceptionDescription exceptionDescription;
          exceptionDescription << "No transportation process found ";
          G4Exception("G4ITStepProcessor::DoDefinePhysicalStepLength",
                      "ITStepProcessor0009",
                      FatalErrorInArgument,
                      exceptionDescription);
        }

        fTimeStep = fpTransportation->GetInteractionTimeLeft();

        if(fpTrack->GetNextVolume() != 0) fpState->fStepStatus = fGeomBoundary;
        else fpState->fStepStatus = fWorldBoundary;
      }
    }
    else
    {
      if(kp == fpProcessInfo->MAXofAlongStepLoops - 1)
      {
        fpTransportation = dynamic_cast<G4ITTransportation*>(fpCurrentProcess);

        if(!fpTransportation)
        {
          G4ExceptionDescription exceptionDescription;
          exceptionDescription << "No transportation process found ";
          G4Exception("G4ITStepProcessor::DoDefinePhysicalStepLength",
                      "ITStepProcessor0010",
                      FatalErrorInArgument,
                      exceptionDescription);
        }

        fTimeStep = fpTransportation->GetInteractionTimeLeft();
      }
    }

    // Handle PostStep processes sending back time steps rather than space length
    if(proposedTimeStep < fTimeStep)
    {
      if(fPostStepAtTimeDoItProcTriggered < fpProcessInfo->MAXofPostStepLoops)
      {
        if((fpState->fSelectedPostStepDoItVector)[fPostStepAtTimeDoItProcTriggered] == InActivated)
        {
          (fpState->fSelectedPostStepDoItVector)[fPostStepAtTimeDoItProcTriggered] =
              NotForced;
          // (fpState->fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] = InActivated;

          fpState->fStepStatus = fPostStepDoItProc;
          fpStep->GetPostStepPoint()->SetProcessDefinedStep(processWithPostStepGivenByTimeStep);

          fTimeStep = proposedTimeStep;

          fpTransportation->ComputeStep(*fpTrack,
                                        *fpStep,
                                        fTimeStep,
                                        fpState->fPhysicalStep);
        }
      }
    }
    else
    {
      if(fPostStepDoItProcTriggered < fpProcessInfo->MAXofPostStepLoops)
      {
        if((fpState->fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] == InActivated)
        {
          (fpState->fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] =
              NotForced;
        }
      }
    }

//  fpCurrentProcess->SetProcessState(0);
    fpCurrentProcess->ResetProcessState();

    // Make sure to check the safety, even if Step is not limited
    //  by this process.                      J. Apostolakis, June 20, 1998
    //
    if(safetyProposedToAndByProcess < fpState->fProposedSafety)
    // proposedSafety keeps the smallest value:
    fpState->fProposedSafety = safetyProposedToAndByProcess;
    else
    // safetyProposedToAndByProcess always proposes a valid safety:
    safetyProposedToAndByProcess = fpState->fProposedSafety;

  }

  fpITrack->GetTrackingInfo()->SetNavigatorState(fpNavigator->GetNavigatorState());
  fpNavigator->ResetNavigatorState();
}

//______________________________________________________________________________
