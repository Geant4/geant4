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
// G4SteppingManager class implementation
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------

#include "G4SteppingManager.hh"
#include "G4SteppingVerbose.hh"
#include "G4SteppingVerboseWithUnits.hh"
#include "G4UImanager.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4SteppingControl.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4VSensitiveDetector.hh"    // Include from 'hits/digi'
#include "G4GeometryTolerance.hh"
#include "G4ParticleTable.hh"
#include "G4Profiler.hh"
#include "G4TiMemory.hh"

// #define debug

//////////////////////////////////////
G4SteppingManager::G4SteppingManager()
//////////////////////////////////////
{
  // Construct simple 'has-a' related objects

  fStep = new G4Step();
  fSecondary = fStep->NewSecondaryVector();
  fPreStepPoint  = fStep->GetPreStepPoint();
  fPostStepPoint = fStep->GetPostStepPoint();

#ifdef G4VERBOSE
  fVerbose = G4VSteppingVerbose::GetInstance();
  if(fVerbose==nullptr)
  {
    if(G4VSteppingVerbose::GetMasterInstance()==nullptr)
    {
      G4int prec = G4SteppingVerbose::BestUnitPrecision();
      if(prec > 0)
      { fVerbose = new G4SteppingVerboseWithUnits(prec); }
      else
      { fVerbose = new G4SteppingVerbose(); }
    }
    else
    { fVerbose = G4VSteppingVerbose::GetMasterInstance()->Clone(); }
    KillVerbose = true;
  }
  else
  { KillVerbose = false; }
  fVerbose -> SetManager(this);
#endif

   SetNavigator(G4TransportationManager::GetTransportationManager()
                ->GetNavigatorForTracking());

   fSelectedAtRestDoItVector
      = new G4SelectedAtRestDoItVector(SizeOfSelectedDoItVector,0);
   fSelectedAlongStepDoItVector
      = new G4SelectedAlongStepDoItVector(SizeOfSelectedDoItVector,0);
   fSelectedPostStepDoItVector
      = new G4SelectedPostStepDoItVector(SizeOfSelectedDoItVector,0);

   SetNavigator(G4TransportationManager::GetTransportationManager()
     ->GetNavigatorForTracking());

   physIntLength = DBL_MAX; 
   kCarTolerance = 0.5*G4GeometryTolerance::GetInstance()
                 ->GetSurfaceTolerance();

   fNoProcess = new G4NoProcess;
}

///////////////////////////////////////
G4SteppingManager::~G4SteppingManager()
///////////////////////////////////////
{
  fTouchableHandle = 0;

  // Destruct simple 'has-a' objects
  //
  fStep->DeleteSecondaryVector();

  // delete fSecondary;
  delete fStep;
  delete fSelectedAtRestDoItVector;
  delete fSelectedAlongStepDoItVector;
  delete fSelectedPostStepDoItVector;
  delete fUserSteppingAction;
  #ifdef G4VERBOSE
    if(KillVerbose) delete fVerbose;
  #endif
}

//////////////////////////////////////////
G4StepStatus G4SteppingManager::Stepping()
//////////////////////////////////////////
{
#ifdef GEANT4_USE_TIMEMORY
  ProfilerConfig profiler{ fStep };
#endif

  //--------
  // Prelude
  //--------
  #ifdef G4VERBOSE
    if(verboseLevel>0)
    {
      fVerbose->NewStep();
    }
    else if (verboseLevel==-1)
    { 
      G4VSteppingVerbose::SetSilent(1);
    }
    else
    {
      G4VSteppingVerbose::SetSilent(0);
    }
  #endif 

  // Store last PostStepPoint to PreStepPoint, and swap current and nex
  // volume information of G4Track. Reset total energy deposit in one Step.
  //
  fStep->CopyPostToPreStepPoint();
  fStep->ResetTotalEnergyDeposit();

  // Switch next touchable in track to current one
  //
  fTrack->SetTouchableHandle(fTrack->GetNextTouchableHandle());

  // Reset the secondary particles
  //
  fN2ndariesAtRestDoIt = 0;
  fN2ndariesAlongStepDoIt = 0;
  fN2ndariesPostStepDoIt = 0;

  // Set the volume before it is used (in DefineStepLength() for User Limit)
  //
  fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();

  // Reset the step's auxiliary points vector pointer
  //
  fStep->SetPointerToVectorOfAuxiliaryPoints(nullptr);

  //-----------------
  // AtRest Processes
  //-----------------

  if( fTrack->GetTrackStatus() == fStopButAlive )
  {
    if( MAXofAtRestLoops>0 )
    {
      InvokeAtRestDoItProcs();
      fStepStatus = fAtRestDoItProc;
      fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );
       
      #ifdef G4VERBOSE
        if(verboseLevel>0) fVerbose->AtRestDoItInvoked();
      #endif 

     }
     // Make sure the track is killed
     //
     fTrack->SetTrackStatus( fStopAndKill );
  }

  //---------------------------------
  // AlongStep and PostStep Processes
  //---------------------------------

  else
  {
    // Find minimum Step length demanded by active disc./cont. processes
    DefinePhysicalStepLength();

    // Store the Step length (geometrical length) to G4Step and G4Track
    fStep->SetStepLength( PhysicalStep );
    fTrack->SetStepLength( PhysicalStep );
    G4double GeomStepLength = PhysicalStep;

    // Store StepStatus to PostStepPoint
    fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );

    // Invoke AlongStepDoIt 
    InvokeAlongStepDoItProcs();

    // Get StepStatus from PostStepPoint - a process such as transportation
    // might have changed it.
    fStepStatus = fStep->GetPostStepPoint()->GetStepStatus();

    // Update track by taking into account all changes by AlongStepDoIt
    fStep->UpdateTrack();

    // Update safety after invocation of all AlongStepDoIts
    endpointSafOrigin= fPostStepPoint->GetPosition();
    // endpointSafety=  std::max( proposedSafety - GeomStepLength, 0.);
    endpointSafety=  std::max( proposedSafety - GeomStepLength, kCarTolerance);

    fStep->GetPostStepPoint()->SetSafety( endpointSafety );

    #ifdef G4VERBOSE
      if(verboseLevel>0) fVerbose->AlongStepDoItAllDone();
    #endif

    // Invoke PostStepDoIt
    InvokePostStepDoItProcs();

    #ifdef G4VERBOSE
      if(verboseLevel>0) fVerbose->PostStepDoItAllDone();
    #endif
  }

  //-------
  // Finale
  //-------

  // Update 'TrackLength' and remeber the Step length of the current Step
  //
  fTrack->AddTrackLength(fStep->GetStepLength());
  fPreviousStepSize = fStep->GetStepLength();
  fStep->SetTrack(fTrack);

  #ifdef G4VERBOSE
    if(verboseLevel>0) fVerbose->StepInfo();
  #endif

  // Send G4Step information to Hit/Dig if the volume is sensitive
  //
  fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
  StepControlFlag =  fStep->GetControlFlag();
  if( fCurrentVolume != nullptr && StepControlFlag != AvoidHitInvocation )
  {
    fSensitive = fStep->GetPreStepPoint()->GetSensitiveDetector();
    if( fSensitive != 0 )
    {
      fSensitive->Hit(fStep);
    }
  }

  // User intervention process
  //
  if( fUserSteppingAction != nullptr )
  {
    fUserSteppingAction->UserSteppingAction(fStep);
  }

  G4UserSteppingAction* regionalAction =
    fCurrentVolume->GetLogicalVolume()->GetRegion()->GetRegionalSteppingAction();

  if(regionalAction)
    regionalAction->UserSteppingAction(fStep);

  // Stepping process finish. Return the value of the StepStatus
  //
  return fStepStatus;
}

///////////////////////////////////////////////////////////
void G4SteppingManager::SetInitialStep(G4Track* valueTrack)
///////////////////////////////////////////////////////////
{
  // Set up several local variables
  //
  PreStepPointIsGeom = false;
  FirstStep = true;
  fParticleChange = nullptr;
  fPreviousStepSize = 0.;
  fStepStatus = fUndefined;

  fTrack = valueTrack;
  Mass = fTrack->GetDynamicParticle()->GetMass();

  PhysicalStep = 0.;
  GeometricalStep = 0.;
  CorrectedStep = 0.;
  PreStepPointIsGeom = false;
  FirstStep = false;

  TempInitVelocity = 0.;
  TempVelocity = 0.;
  sumEnergyChange = 0.;

  // If the primary track has 'Suspend' or 'PostponeToNextEvent' state,
  // set the track state to 'Alive'
  //
  if ( (fTrack->GetTrackStatus()==fSuspend)
    || (fTrack->GetTrackStatus()==fPostponeToNextEvent) )
  { 
    fTrack->SetTrackStatus(fAlive);
  }

  // If the primary track has 'zero' kinetic energy, set the track
  // state to 'StopButAlive'
  //
  if(fTrack->GetKineticEnergy() <= 0.0)
  {
    fTrack->SetTrackStatus( fStopButAlive );
  }

  // Set Touchable to track and a private attribute of G4SteppingManager
 
  if ( ! fTrack->GetTouchableHandle() )
  {
    G4ThreeVector direction= fTrack->GetMomentumDirection();
    fNavigator->LocateGlobalPointAndSetup( fTrack->GetPosition(),
                                           &direction, false, false );
    fTouchableHandle = fNavigator->CreateTouchableHistory();
    fTrack->SetTouchableHandle( fTouchableHandle );
    fTrack->SetNextTouchableHandle( fTouchableHandle );
  }
  else
  {
    fTrack->SetNextTouchableHandle( fTouchableHandle = fTrack->GetTouchableHandle() );
    G4VPhysicalVolume* oldTopVolume = fTrack->GetTouchableHandle()->GetVolume();
    G4VPhysicalVolume* newTopVolume =
      fNavigator->ResetHierarchyAndLocate( fTrack->GetPosition(), 
                     fTrack->GetMomentumDirection(),
                     *((G4TouchableHistory*)fTrack->GetTouchableHandle()()) );
    if ( newTopVolume != oldTopVolume
      || oldTopVolume->GetRegularStructureId() == 1 )
    { 
      fTouchableHandle = fNavigator->CreateTouchableHistory();
      fTrack->SetTouchableHandle( fTouchableHandle );
      fTrack->SetNextTouchableHandle( fTouchableHandle );
    }
  }

  // Set OriginTouchableHandle for primary track
  //
  if(fTrack->GetParentID()==0)
  {
    fTrack->SetOriginTouchableHandle(fTrack->GetTouchableHandle());
  }

  // Set vertex information of G4Track at here
  //
  if ( fTrack->GetCurrentStepNumber() == 0 )
  {
    fTrack->SetVertexPosition( fTrack->GetPosition() );
    fTrack->SetVertexMomentumDirection( fTrack->GetMomentumDirection() );
    fTrack->SetVertexKineticEnergy( fTrack->GetKineticEnergy() );
    fTrack->SetLogicalVolumeAtVertex( fTrack->GetVolume()->GetLogicalVolume() );
  }

  // Initial set up for attributes of 'G4SteppingManager'
  fCurrentVolume = fTouchableHandle->GetVolume();

  // If track is already outside the world boundary, kill it
  //
  if( fCurrentVolume==nullptr )
  {
    // If the track is a primary, stop processing
    if(fTrack->GetParentID()==0)
    {
      G4cerr << "ERROR - G4SteppingManager::SetInitialStep()" << G4endl
             << "        Primary particle starting at - "
             << fTrack->GetPosition()
             << " - is outside of the world volume." << G4endl;
      G4Exception("G4SteppingManager::SetInitialStep()", "Tracking0010",
                  FatalException, "Primary vertex outside of the world!");
    }

    fTrack->SetTrackStatus( fStopAndKill );
    G4cout << "WARNING - G4SteppingManager::SetInitialStep()" << G4endl
           << "          Initial track position is outside world! - "
           << fTrack->GetPosition() << G4endl;
   }
   else
   {
     // Initial set up for attributes of 'Step'
     fStep->InitializeStep( fTrack );
   }

   #ifdef G4VERBOSE
     if(verboseLevel>0) fVerbose->TrackingStarted();
   #endif
}

/////////////////////////////////////////////////
void G4SteppingManager::GetProcessNumber()
/////////////////////////////////////////////////
{
  #ifdef debug
    G4cout << "G4SteppingManager::GetProcessNumber: is called track="
           << fTrack << G4endl;
  #endif

  G4ProcessManager* pm= fTrack->GetDefinition()->GetProcessManager();
  if(pm == nullptr)
  {
    G4cerr << "ERROR - G4SteppingManager::GetProcessNumber()" << G4endl
           << "        ProcessManager is NULL for particle = "
           << fTrack->GetDefinition()->GetParticleName() << ", PDG_code = "
           << fTrack->GetDefinition()->GetPDGEncoding() << G4endl;
    G4Exception("G4SteppingManager::GetProcessNumber()", "Tracking0011",
                FatalException, "Process Manager is not found.");
    return;
  }

  // AtRestDoits
  //
  MAXofAtRestLoops =        pm->GetAtRestProcessVector()->entries();
  fAtRestDoItVector =       pm->GetAtRestProcessVector(typeDoIt);
  fAtRestGetPhysIntVector = pm->GetAtRestProcessVector(typeGPIL);

  #ifdef debug
    G4cout << "G4SteppingManager::GetProcessNumber: #ofAtRest="
           << MAXofAtRestLoops << G4endl;
  #endif

  // AlongStepDoits
  //
  MAXofAlongStepLoops = pm->GetAlongStepProcessVector()->entries();
  fAlongStepDoItVector = pm->GetAlongStepProcessVector(typeDoIt);
  fAlongStepGetPhysIntVector = pm->GetAlongStepProcessVector(typeGPIL);

  #ifdef debug
    G4cout << "G4SteppingManager::GetProcessNumber:#ofAlongStp="
           << MAXofAlongStepLoops << G4endl;
  #endif

  // PostStepDoits
  //
  MAXofPostStepLoops = pm->GetPostStepProcessVector()->entries();
  fPostStepDoItVector = pm->GetPostStepProcessVector(typeDoIt);
  fPostStepGetPhysIntVector = pm->GetPostStepProcessVector(typeGPIL);

  #ifdef debug
    G4cout << "G4SteppingManager::GetProcessNumber: #ofPostStep="
           << MAXofPostStepLoops << G4endl;
  #endif

  if (SizeOfSelectedDoItVector<MAXofAtRestLoops    ||
      SizeOfSelectedDoItVector<MAXofAlongStepLoops ||
      SizeOfSelectedDoItVector<MAXofPostStepLoops  )
  {
    G4cerr << "ERROR - G4SteppingManager::GetProcessNumber()" << G4endl
           << "        SizeOfSelectedDoItVector= " << SizeOfSelectedDoItVector
           << " ; is smaller then one of MAXofAtRestLoops= "
           << MAXofAtRestLoops << G4endl
           << "        or MAXofAlongStepLoops= " << MAXofAlongStepLoops
           << " or MAXofPostStepLoops= " << MAXofPostStepLoops << G4endl;
    G4Exception("G4SteppingManager::GetProcessNumber()",
                "Tracking0012", FatalException,
                "The array size is smaller than the actual No of processes.");
  }
}

// ************************************************************************
//
//  Private Member Functions
//
// ************************************************************************

/////////////////////////////////////////////////////////
void G4SteppingManager::DefinePhysicalStepLength()
/////////////////////////////////////////////////////////
{
  // ReSet the counter etc.
  //
  PhysicalStep  = DBL_MAX;          // Initialize by a huge number    
  physIntLength = DBL_MAX;          // Initialize by a huge number    

  #ifdef G4VERBOSE
    if(verboseLevel>0) fVerbose->DPSLStarted();
  #endif

  // GPIL for PostStep
  //
  fPostStepDoItProcTriggered = MAXofPostStepLoops;

  for(std::size_t np=0; np<MAXofPostStepLoops; ++np)
  {
    fCurrentProcess = (*fPostStepGetPhysIntVector)((G4int)np);
    if (fCurrentProcess == nullptr)
    {
      (*fSelectedPostStepDoItVector)[np] = InActivated;
      continue;
    } // NULL means the process is inactivated by a user on fly

    physIntLength = fCurrentProcess->PostStepGPIL( *fTrack,
                                     fPreviousStepSize, &fCondition );
    #ifdef G4VERBOSE
      if(verboseLevel>0) fVerbose->DPSLPostStep();
    #endif

    switch (fCondition)
    {
      case ExclusivelyForced:
        (*fSelectedPostStepDoItVector)[np] = ExclusivelyForced;
        fStepStatus = fExclusivelyForcedProc;
        fStep->GetPostStepPoint()->SetProcessDefinedStep(fCurrentProcess);
        break;
      case Conditionally:
        // (*fSelectedPostStepDoItVector)[np] = Conditionally;
        G4Exception("G4SteppingManager::DefinePhysicalStepLength()",
                    "Tracking1001", FatalException,
                    "This feature no more supported");
        break;
      case Forced:
        (*fSelectedPostStepDoItVector)[np] = Forced;
        break;
      case StronglyForced:
        (*fSelectedPostStepDoItVector)[np] = StronglyForced;
        break;
      default:
        (*fSelectedPostStepDoItVector)[np] = InActivated;
        break;
    }

    if (fCondition==ExclusivelyForced)
    { 
      for(std::size_t nrest=np+1; nrest<MAXofPostStepLoops; ++nrest)
      { 
        (*fSelectedPostStepDoItVector)[nrest] = InActivated; 
      } 
      return;  // Take note the 'return' at here !!! 
    } 
    else
    {
      if(physIntLength < PhysicalStep )
      {
        PhysicalStep = physIntLength;
        fStepStatus = fPostStepDoItProc;
        fPostStepDoItProcTriggered = G4int(np);
        fStep->GetPostStepPoint()->SetProcessDefinedStep(fCurrentProcess);
      }
    }
  }

  if (fPostStepDoItProcTriggered<MAXofPostStepLoops)
  {
    if ((*fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] == InActivated)
    {
      (*fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] = NotForced;
    }
  }

  // GPIL for AlongStep
  //
  proposedSafety = DBL_MAX;
  G4double safetyProposedToAndByProcess = proposedSafety;
  G4bool delegateToTransportation = false;

  for(std::size_t kp=0; kp<MAXofAlongStepLoops; ++kp)
  {
    fCurrentProcess = (*fAlongStepGetPhysIntVector)[(G4int)kp];
    if (fCurrentProcess == nullptr) continue;
      // NULL means the process is inactivated by a user on fly

    physIntLength = fCurrentProcess->AlongStepGPIL( *fTrack,
                                     fPreviousStepSize, PhysicalStep,
                                     safetyProposedToAndByProcess,
                                     &fGPILSelection );
    #ifdef G4VERBOSE
      if(verboseLevel>0) fVerbose->DPSLAlongStep();
    #endif

    if(physIntLength < PhysicalStep)
    {
      PhysicalStep = physIntLength;

      // Check if the process wants to be the GPIL winner. For example,
      // multi-scattering proposes Step limit, but won't be the winner
      //
      if(fGPILSelection==CandidateForSelection)
      {
        fStepStatus = fAlongStepDoItProc;
        fStep->GetPostStepPoint()->SetProcessDefinedStep(fCurrentProcess);
      }
      else if(fCurrentProcess->GetProcessType()==fParallel)
      { // a parallel world is proposing the shortest but expecting Transportation
        // to win.
        delegateToTransportation = true;
      }

      // Transportation is assumed to be the last process in the vector
      // Transportation is winning
      if(kp == MAXofAlongStepLoops-1) 
      {
        // This used to set fStepStatus = fGeomBoundary, but it was moved to
        // G4Transportation::AlongStepDoIt where the process can actually
        // decide if there is a volume boundary.
        delegateToTransportation = false;
      }
    }

    // Make sure to check the safety, even if Step is not limited 
    // by this process
    // 
    if (safetyProposedToAndByProcess < proposedSafety)
    {
      // proposedSafety keeps the smallest value
      //
      proposedSafety = safetyProposedToAndByProcess;
    }
    else
    {
      // safetyProposedToAndByProcess always proposes a valid safety
      //
      safetyProposedToAndByProcess = proposedSafety;
    }
  } 
  if(delegateToTransportation)
  {
    fStepStatus = fGeomBoundary;
    fStep->GetPostStepPoint()->SetProcessDefinedStep(fCurrentProcess);
  }
}

//////////////////////////////////////////////////////
G4int G4SteppingManager::ProcessSecondariesFromParticleChange()
//////////////////////////////////////////////////////
{
  G4Track* tempSecondaryTrack;
  G4int    num2ndaries;
  G4int    pushedSecondaries = 0;

  num2ndaries = fParticleChange->GetNumberOfSecondaries();
  if(num2ndaries == 0)
  {
    return 0;
  }

  // Get the creator process. This may be different from fCurrentProcess for a
  // "combined" process such as G4GammaGeneralProcess.
  const G4VProcess* creatorProcess = fCurrentProcess->GetCreatorProcess();

  for(G4int DSecLoop=0; DSecLoop< num2ndaries; ++DSecLoop)
  {
    tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);

    // Set parentID
    tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );

    // Set the process pointer which created this track
    tempSecondaryTrack->SetCreatorProcess( creatorProcess );

    // If this 2ndry particle has 'zero' kinetic energy, make sure
    // it invokes a rest process at the beginning of the tracking
    //
    if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN)
    {
      G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()
                             ->GetProcessManager();
      if(pm == nullptr)
      {
        G4ExceptionDescription ED;
        ED << "A track without proper process manager is pushed\n"
           << "into the track stack.\n"
           << " Particle name : "
           << tempSecondaryTrack->GetDefinition()->GetParticleName()
           << " -- created by " << creatorProcess->GetProcessName() << ".";
        G4Exception("G4SteppingManager::ProcessSecondariesFromParticleChange()",
                    "Tracking10051", FatalException, ED);
      }
      if (pm->GetAtRestProcessVector()->entries()>0)
      {
        tempSecondaryTrack->SetTrackStatus( fStopButAlive );
        fSecondary->push_back( tempSecondaryTrack );
        ++pushedSecondaries;
      }
      else
      {
        delete tempSecondaryTrack;
      }
    }
    else
    {
      fSecondary->push_back( tempSecondaryTrack );
      ++pushedSecondaries;
    }
  } //end of loop on secondary

  return pushedSecondaries;
}

//////////////////////////////////////////////////////
void G4SteppingManager::InvokeAtRestDoItProcs()
//////////////////////////////////////////////////////
{
  // Select the rest process which has the shortest time before
  // it is invoked. In rest processes, GPIL()
  // returns the time before a process occurs

  G4double lifeTime, shortestLifeTime;

  fAtRestDoItProcTriggered = 0;
  shortestLifeTime = DBL_MAX;

  for( std::size_t ri=0 ; ri < MAXofAtRestLoops ; ++ri )
  {
    fCurrentProcess = (*fAtRestGetPhysIntVector)[(G4int)ri];
    if (fCurrentProcess == nullptr)
    {
      (*fSelectedAtRestDoItVector)[ri] = InActivated;
      continue;
    }   // nullptr means the process is inactivated by a user on fly

    lifeTime = fCurrentProcess->AtRestGPIL( *fTrack, &fCondition );

    if(fCondition == Forced)
    {
      (*fSelectedAtRestDoItVector)[ri] = Forced;
    }
    else
    {
      (*fSelectedAtRestDoItVector)[ri] = InActivated;
      if(lifeTime < shortestLifeTime )
      {
        shortestLifeTime = lifeTime;
        fAtRestDoItProcTriggered = G4int(ri);
        fStep->GetPostStepPoint()->SetProcessDefinedStep(fCurrentProcess);
      }
    }
  }

  (*fSelectedAtRestDoItVector)[fAtRestDoItProcTriggered] = NotForced;

  fStep->SetStepLength( 0. );  // the particle has stopped
  fTrack->SetStepLength( 0. );

  // Condition to avoid that stable ions are handled by Radioactive Decay.
  // We use a very large time threshold (many orders of magnitude bigger than
  // the universe's age) but not DBL_MAX because shortestLifeTime can be
  // sometimes slightly smaller for stable ions.
  if(shortestLifeTime < 1.0e+100)  // Unstable ion at rest: Radioactive Decay will decay it
  {
    // invoke selected process
    //
    for(std::size_t np=0; np<MAXofAtRestLoops; ++np)
    {
      //
      // Note: DoItVector has inverse order against GetPhysIntVector
      //       and SelectedAtRestDoItVector.
      //
      if( (*fSelectedAtRestDoItVector)[MAXofAtRestLoops-np-1] != InActivated)
      {
        fCurrentProcess = (*fAtRestDoItVector)[(G4int)np];
        fParticleChange = fCurrentProcess->AtRestDoIt(*fTrack, *fStep);

        // Update Step
        //
        fParticleChange->UpdateStepForAtRest(fStep);

        // Now Store the secondaries from ParticleChange to SecondaryList
        fN2ndariesAtRestDoIt += ProcessSecondariesFromParticleChange();

        // clear ParticleChange
        fParticleChange->Clear();

      }  // if(fSelectedAtRestDoItVector[np] != InActivated){
    }  // for(std::size_t np=0; np<MAXofAtRestLoops; ++np){
  }
  else  // Stable ion at rest
  {
    fStep->GetPostStepPoint()->SetProcessDefinedStep( fNoProcess );
  }  // if(shortestLifeTime < 1.0e+100)

  fStep->UpdateTrack();

  fTrack->SetTrackStatus( fStopAndKill );
}

/////////////////////////////////////////////////////////
void G4SteppingManager::InvokeAlongStepDoItProcs()
/////////////////////////////////////////////////////////
{
  // If the current Step is defined by a 'ExclusivelyForced' 
  // PostStepDoIt, then don't invoke any AlongStepDoIt
  //
  if(fStepStatus == fExclusivelyForcedProc)
  {
    return;               // Take note 'return' is here !!!
  }

  // Invoke all active continuous processes
  //
  for( std::size_t ci=0; ci<MAXofAlongStepLoops; ++ci )
  {
    fCurrentProcess = (*fAlongStepDoItVector)[(G4int)ci];
    if (fCurrentProcess== 0) continue;
      // NULL means the process is inactivated by a user on fly.

    fParticleChange = fCurrentProcess->AlongStepDoIt( *fTrack, *fStep );

    // Update the PostStepPoint of Step according to ParticleChange
    fParticleChange->UpdateStepForAlongStep(fStep);

    #ifdef G4VERBOSE
      if(verboseLevel>0) fVerbose->AlongStepDoItOneByOne();
    #endif

    // Now Store the secondaries from ParticleChange to SecondaryList
    fN2ndariesAlongStepDoIt += ProcessSecondariesFromParticleChange();
     
    // Set the track status according to what the process defined
    // if kinetic energy >0, otherwise set  fStopButAlive
    //
    fTrack->SetTrackStatus( fParticleChange->GetTrackStatus() );
     
    // clear ParticleChange
    fParticleChange->Clear();
  }

  fStep->UpdateTrack();
  G4TrackStatus fNewStatus = fTrack->GetTrackStatus();

  if ( fNewStatus == fAlive && fTrack->GetKineticEnergy() <= DBL_MIN )
  {
    if(MAXofAtRestLoops>0) fNewStatus = fStopButAlive;
    else                   fNewStatus = fStopAndKill;
    fTrack->SetTrackStatus( fNewStatus );
  }
}

////////////////////////////////////////////////////////
void G4SteppingManager::InvokePostStepDoItProcs()
////////////////////////////////////////////////////////
{
  // Invoke the specified discrete processes
  //
  for(std::size_t np=0; np<MAXofPostStepLoops; ++np)
  {
    //
    // Note: DoItVector has inverse order against GetPhysIntVector
    //       and SelectedPostStepDoItVector.
    //
    G4int Cond = (*fSelectedPostStepDoItVector)[MAXofPostStepLoops-np-1];
    if(Cond != InActivated)
    {
      if( ((Cond == NotForced) && (fStepStatus == fPostStepDoItProc)) ||
          ((Cond == Forced) && (fStepStatus != fExclusivelyForcedProc)) ||
          ((Cond == ExclusivelyForced) && (fStepStatus == fExclusivelyForcedProc)) || 
          ((Cond == StronglyForced) ) )
      {
        InvokePSDIP(np);
        if ((np==0) && (fTrack->GetNextVolume() == nullptr))
        {
          fStepStatus = fWorldBoundary;
          fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );
        }
      }
    }

    // Exit from PostStepLoop if the track has been killed,
    // but extra treatment for processes with Strongly Forced flag
    //
    if(fTrack->GetTrackStatus() == fStopAndKill)
    {
      for(std::size_t np1=np+1; np1<MAXofPostStepLoops; ++np1)
      { 
        G4int Cond2 = (*fSelectedPostStepDoItVector)[MAXofPostStepLoops-np1-1];
        if (Cond2 == StronglyForced)
        {
          InvokePSDIP(np1);
        }
      }
      break;
    }
  }
}

////////////////////////////////////////////////////////
void G4SteppingManager::InvokePSDIP(size_t np)
////////////////////////////////////////////////////////
{
  fCurrentProcess = (*fPostStepDoItVector)[(G4int)np];
  fParticleChange = fCurrentProcess->PostStepDoIt( *fTrack, *fStep);

  // Update PostStepPoint of Step according to ParticleChange
  fParticleChange->UpdateStepForPostStep(fStep);

  #ifdef G4VERBOSE
    if(verboseLevel>0) fVerbose->PostStepDoItOneByOne();
  #endif

  // Update G4Track according to ParticleChange after each PostStepDoIt
  fStep->UpdateTrack();

  // Update safety after each invocation of PostStepDoIts
  fStep->GetPostStepPoint()->SetSafety( CalculateSafety() );

  // Now Store the secondaries from ParticleChange to SecondaryList
  fN2ndariesPostStepDoIt += ProcessSecondariesFromParticleChange();

  // Set the track status according to what the process defined
  fTrack->SetTrackStatus( fParticleChange->GetTrackStatus() );

  // clear ParticleChange
  fParticleChange->Clear();
}
