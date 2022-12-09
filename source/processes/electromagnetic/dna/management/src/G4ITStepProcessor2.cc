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
#include "G4LossTableManager.hh"
#include "G4EnergyLossTables.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VITProcess.hh"
#include "G4TrackingInformation.hh"
#include "G4IT.hh"
#include "G4ITTrackingManager.hh"
#include "G4ITTransportation.hh"

#include "G4ITNavigator.hh"             // Include from 'geometry'

#include "G4ITSteppingVerbose.hh"
#include "G4VITSteppingVerbose.hh"

#include "G4ITTrackHolder.hh"
#include "G4ITReaction.hh"

#include "G4TrackingInformation.hh"

//#define DEBUG_MEM 1

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
using G4MemStat::MemStat;
#endif

void G4ITStepProcessor::DealWithSecondaries(G4int& counter)
{
  // Now Store the secondaries from ParticleChange to SecondaryList
  G4Track* tempSecondaryTrack;

  for(G4int DSecLoop = 0; DSecLoop < fpParticleChange->GetNumberOfSecondaries();
      DSecLoop++)
  {
    tempSecondaryTrack = fpParticleChange->GetSecondary(DSecLoop);

    if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
    {
      ApplyProductionCut(tempSecondaryTrack);
    }

    // Set parentID
    tempSecondaryTrack->SetParentID(fpTrack->GetTrackID());

    // Set the process pointer which created this track
    tempSecondaryTrack->SetCreatorProcess(fpCurrentProcess);

    // If this 2ndry particle has 'zero' kinetic energy, make sure
    // it invokes a rest process at the beginning of the tracking
    if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN)
    {
      G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
      if (pm->GetAtRestProcessVector()->entries()>0)
      {
        tempSecondaryTrack->SetTrackStatus( fStopButAlive );
        fpSecondary->push_back( tempSecondaryTrack );
        fN2ndariesAtRestDoIt++;
      }
      else
      {
        delete tempSecondaryTrack;
      }
    }
    else
    {
      fpSecondary->push_back( tempSecondaryTrack );
      counter++;
    }
  } //end of loop on secondary
}

//_________________________________________________________________________

void G4ITStepProcessor::DoIt(double timeStep)

// Call the process having the min step length or just propagate the track on the given time step

// If the track is "leading the step" (ie one of its process has been selected
// as the one having the minimum time step over all tracks and processes),
// it will undergo its selected processes. Otherwise, it will just propagate the track
// on the given time step.

{
  if(fpVerbose) fpVerbose->DoItStarted();

  G4TrackManyList* mainList = fpTrackContainer->GetMainList();
  G4TrackManyList::iterator it = mainList->end();
  it--;
  std::size_t initialSize = mainList->size();

//    G4cout << "initialSize = " << initialSize << G4endl;

  for(std::size_t i = 0 ; i < initialSize ; ++i)
  {

//      G4cout << "i = " << i << G4endl;

    G4Track* track = *it;
    if (!track)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No track was pop back the main track list.";
      G4Exception("G4ITStepProcessor::DoIt", "NO_TRACK",
                  FatalException, exceptionDescription);
    }
    // G4TrackManyList::iterator next_it (it);
    // next_it--;
    // it = next_it;

    it--;
    // Must be called before EndTracking(track)
    // Otherwise the iterator will point to the list of killed tracks

    if(track->GetTrackStatus() == fStopAndKill)
    {
      fpTrackingManager->EndTracking(track);
//      G4cout << GetIT(track)->GetName() << G4endl;
//      G4cout << " ************************ CONTINUE ********************" << G4endl;
      continue;
    }

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    MemStat mem_first, mem_second, mem_diff;
    mem_first = MemoryUsage();
#endif

    Stepping(track, timeStep);

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    MemStat mem_intermediaire = MemoryUsage();
    mem_diff = mem_intermediaire-mem_first;
    G4cout << "\t\t >> || MEM || In DoIT with track "
        << track->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    ExtractDoItData();

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    mem_second = MemoryUsage();
    mem_diff = mem_second-mem_first;
    G4cout << "\t >> || MEM || In DoIT with track "
        << track->GetTrackID()
        << ", diff is : " << mem_diff << G4endl;
#endif
  }


  fpTrackContainer->MergeSecondariesWithMainList();
  fpTrackContainer->KillTracks(); // (18-06-15 : MK) Remove it ?
  fLeadingTracks.Reset();
}

//_________________________________________________________________________

void G4ITStepProcessor::ExtractDoItData()
{
  if (!fpTrack)
  {
    CleanProcessor();
    return;
  }

  G4TrackStatus status = fpTrack->GetTrackStatus();

  switch (status)
  {
    case fAlive:
    case fStopButAlive:
    case fSuspend:
    case fPostponeToNextEvent:
    default:
      PushSecondaries();
      break;

    case fStopAndKill:
      G4ITReactionSet::Instance()->RemoveReactionSet(fpTrack);
      PushSecondaries();
//      G4TrackList::Pop(fpTrack);
      fpTrackingManager->EndTracking(fpTrack);
//      fTrackContainer->PushToKill(fpTrack);
      break;

    case fKillTrackAndSecondaries:
      G4ITReactionSet::Instance()->RemoveReactionSet(fpTrack);
      if (fpSecondary)
      {
        for (std::size_t i = 0; i < fpSecondary->size(); ++i)
        {
          delete (*fpSecondary)[i];
        }
        fpSecondary->clear();
      }
//      G4TrackList::Pop(fpTrack);
      fpTrackingManager->EndTracking(fpTrack);
//      fTrackContainer->PushToKill(fpTrack);
      break;
  }

  CleanProcessor();
}

//_________________________________________________________________________

void G4ITStepProcessor::PushSecondaries()
{
  if (!fpSecondary || fpSecondary->empty())
  {
    // DEBUG
    //      G4cout << "NO SECONDARIES !!! " << G4endl;
    return;
  }

  // DEBUG
  //    G4cout << "There are secondaries : "<< secondaries -> size() << G4endl ;

  G4TrackVector::iterator secondaries_i = fpSecondary->begin();

  for (; secondaries_i != fpSecondary->end(); ++secondaries_i)
  {
    G4Track* secondary = *secondaries_i;
    fpTrackContainer->_PushTrack(secondary);
  }
}

//______________________________________________________________________________

void G4ITStepProcessor::Stepping(G4Track* track, const double & timeStep)
{

#ifdef DEBUG_MEM
  MemStat mem_first, mem_second, mem_diff;
#endif

#ifdef DEBUG_MEM
  mem_first = MemoryUsage();
#endif

  CleanProcessor();

#ifdef DEBUG_MEM
  MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || After CleanProcessor " << track->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

  if(track == 0) return; // maybe put an exception here
  fTimeStep = timeStep;
  SetTrack(track);
  DoStepping();
}
//______________________________________________________________________________

// ************************************************************************
//  Stepping
// ************************************************************************
void G4ITStepProcessor::DoStepping()
{
  SetupMembers();

#ifdef DEBUG_MEM
  MemStat mem_first, mem_second, mem_diff;
#endif

#ifdef DEBUG_MEM
  mem_first = MemoryUsage();
#endif

#ifdef G4VERBOSE
    if(fpVerbose) fpVerbose->PreStepVerbose(fpTrack);
#endif

  if(!fpProcessInfo)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No process info found for particle :"
                         << fpTrack->GetDefinition()->GetParticleName();
    G4Exception("G4ITStepProcessor::DoStepping",
                "ITStepProcessor0012",
                FatalErrorInArgument,
                exceptionDescription);
    return;
  }
//  else if(fpTrack->GetTrackStatus() == fStopAndKill)
//  {
//    fpState->fStepStatus = fUndefined;
//    return;
//  }

  if(fpProcessInfo->MAXofPostStepLoops == 0 &&
      fpProcessInfo->MAXofAlongStepLoops == 0
     && fpProcessInfo->MAXofAtRestLoops == 0)
  {/*
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No process was found for particle :"
                         << fpTrack->GetDefinition()->GetParticleName();
    G4Exception("G4ITStepProcessor::DoStepping",
                "ITStepProcessorNoProcess",
                JustWarning,
                exceptionDescription);

    fpTrack->SetTrackStatus(fStopAndKill);
    fpState->fStepStatus = fUndefined;*/
    return;
  }

  //--------
  // Prelude
  //--------
#ifdef G4VERBOSE
  // !!!!! Verbose
  if(fpVerbose) fpVerbose->NewStep();
#endif

  //---------------------------------
  // AtRestStep, AlongStep and PostStep Processes
  //---------------------------------

  fpNavigator->SetNavigatorState(fpITrack->GetTrackingInfo()->GetNavigatorState());
//        fpNavigator->ResetHierarchyAndLocate( fpTrack->GetPosition(),
//                                              fpTrack->GetMomentumDirection(),
//                                              *((G4TouchableHistory*)fpTrack->GetTouchableHandle()()) );
//        fpNavigator->SetNavigatorState(fpITrack->GetTrackingInfo()->GetNavigatorState());
  // We reset the navigator state before checking for AtRest
  // in case a AtRest processe would use a navigator info

#ifdef DEBUG_MEM
  MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || G4ITStepProcessor::DoStepping || After dealing with navigator with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

  if(fpTrack->GetTrackStatus() == fStopButAlive)
  {
    if(fpProcessInfo->MAXofAtRestLoops > 0 && fpProcessInfo->fpAtRestDoItVector
        != 0) // second condition to make coverity happy
    {
      //-----------------
      // AtRestStepDoIt
      //-----------------
      InvokeAtRestDoItProcs();
      fpState->fStepStatus = fAtRestDoItProc;
      fpStep->GetPostStepPoint()->SetStepStatus(fpState->fStepStatus);

#ifdef G4VERBOSE
            // !!!!! Verbose
      if(fpVerbose) fpVerbose->AtRestDoItInvoked();
#endif

    }
    // Make sure the track is killed
    // fpTrack->SetTrackStatus(fStopAndKill);
  }
  else // if(fTimeStep > 0.) // Bye, because PostStepIL can return 0 => time =0
  {
    if(fpITrack == 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << " !!! TrackID : " << fpTrack->GetTrackID()
                           << G4endl<< " !!! Track status : "<< fpTrack->GetTrackStatus() << G4endl
      << " !!! Particle Name : "<< fpTrack -> GetDefinition() -> GetParticleName() << G4endl
      << "No G4ITStepProcessor::fpITrack found" << G4endl;

      G4Exception("G4ITStepProcessor::DoStepping",
                  "ITStepProcessor0013",
                  FatalErrorInArgument,
                  exceptionDescription);
      return; // to make coverity happy
    }

    if(fpITrack->GetTrackingInfo()->IsLeadingStep() == false)
    {
      // In case the track has NOT the minimum step length
      // Given the final step time, the transportation
      // will compute the final position of the particle
      fpState->fStepStatus = fPostStepDoItProc;
      fpStep->GetPostStepPoint()->SetProcessDefinedStep(fpTransportation);
      FindTransportationStep();
    }

#ifdef DEBUG_MEM
    mem_intermediaire = MemoryUsage();
    mem_diff = mem_intermediaire-mem_first;
    G4cout << "\t\t\t >> || MEM || G4ITStepProcessor::DoStepping || After FindTransportationStep() with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    // Store the Step length (geometrical length) to G4Step and G4Track
    fpTrack->SetStepLength(fpState->fPhysicalStep);
    fpStep->SetStepLength(fpState->fPhysicalStep);

    G4double GeomStepLength = fpState->fPhysicalStep;

    // Store StepStatus to PostStepPoint
    fpStep->GetPostStepPoint()->SetStepStatus(fpState->fStepStatus);

    // Invoke AlongStepDoIt
    InvokeAlongStepDoItProcs();

#ifdef DEBUG_MEM
    mem_intermediaire = MemoryUsage();
    mem_diff = mem_intermediaire-mem_first;
    G4cout << "\t\t\t >> || MEM || G4ITStepProcessor::DoStepping || After InvokeAlongStepDoItProcs() with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

#ifdef G4VERBOSE
    // !!!!! Verbose
    if(fpVerbose) fpVerbose->AlongStepDoItAllDone();
#endif

    // Update track by taking into account all changes by AlongStepDoIt
    // fpStep->UpdateTrack(); // done in InvokeAlongStepDoItProcs

    // Update safety after invocation of all AlongStepDoIts
    fpState->fEndpointSafOrigin = fpPostStepPoint->GetPosition();

    fpState->fEndpointSafety =
        std::max(fpState->fProposedSafety - GeomStepLength, kCarTolerance);

    fpStep->GetPostStepPoint()->SetSafety(fpState->fEndpointSafety);

    if(GetIT(fpTrack)->GetTrackingInfo()->IsLeadingStep())
    {
      // Invoke PostStepDoIt including G4ITTransportation::PSDI
      InvokePostStepDoItProcs();

#ifdef DEBUG_MEM
      mem_intermediaire = MemoryUsage();
      mem_diff = mem_intermediaire-mem_first;
      G4cout << "\t\t\t >> || MEM || G4ITStepProcessor::DoStepping || After InvokePostStepDoItProcs() with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif
#ifdef G4VERBOSE
    // !!!!! Verbose
    if(fpVerbose) fpVerbose->StepInfoForLeadingTrack();
#endif
    }
    else
    {
      // Only invoke transportation and all other forced processes
      InvokeTransportationProc();
      fpStep->GetPostStepPoint()->SetProcessDefinedStep(fpTransportation);

#ifdef DEBUG_MEM
      mem_intermediaire = MemoryUsage();
      mem_diff = mem_intermediaire-mem_first;
      G4cout << "\t\t\t >> || MEM || G4ITStepProcessor::DoStepping || After InvokeTransportationProc() with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif
    }

#ifdef G4VERBOSE
    // !!!!! Verbose
    if(fpVerbose) fpVerbose->PostStepDoItAllDone();
#endif
  }

  fpNavigator->ResetNavigatorState();

#ifdef DEBUG_MEM
  mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || G4ITStepProcessor::DoStepping || After fpNavigator->SetNavigatorState with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

  //-------
  // Finale
  //-------

  // Update 'TrackLength' and remeber the Step length of the current Step
  fpTrack->AddTrackLength(fpStep->GetStepLength());
  fpTrack->IncrementCurrentStepNumber();

//#ifdef G4VERBOSE
// //    !!!!! Verbose
//  if(fpVerbose) fpVerbose->StepInfo();
//#endif

#ifdef G4VERBOSE
    if(fpVerbose) fpVerbose->PostStepVerbose(fpTrack);
#endif

//    G4cout << " G4ITStepProcessor::DoStepping  -- " <<fpTrack->GetTrackID() << " tps = " << fpTrack->GetGlobalTime() << G4endl;

  // Send G4Step information to Hit/Dig if the volume is sensitive
  /***
   fpCurrentVolume = fpStep->GetPreStepPoint()->GetPhysicalVolume();
   StepControlFlag =  fpStep->GetControlFlag();

   if( fpCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation)
   {
   fpSensitive = fpStep->GetPreStepPoint()->
   GetSensitiveDetector();
   if( fpSensitive != 0 )
   {
   fpSensitive->Hit(fpStep);
   }
   }

   User intervention process.
   if( fpUserSteppingAction != 0 )
   {
   fpUserSteppingAction->UserSteppingAction(fpStep);
   }
   G4UserSteppingAction* regionalAction
   = fpStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetRegion()
   ->GetRegionalSteppingAction();
   if( regionalAction ) regionalAction->UserSteppingAction(fpStep);
   ***/
  fpTrackingManager->AppendStep(fpTrack, fpStep);
  // Stepping process finish. Return the value of the StepStatus.

#ifdef DEBUG_MEM
  MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || End of DoStepping() with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

  //    return fpState->fStepStatus;
}

//______________________________________________________________________________

// ************************************************************************
//  AtRestDoIt
// ************************************************************************

void G4ITStepProcessor::InvokeAtRestDoItProcs()
{
  fpStep->SetStepLength(0.);  //the particle has stopped
  fpTrack->SetStepLength(0.);

  G4SelectedAtRestDoItVector& selectedAtRestDoItVector =
      fpState->fSelectedAtRestDoItVector;

  // invoke selected process
  for(std::size_t np = 0; np < fpProcessInfo->MAXofAtRestLoops; ++np)
  {
    //
    // Note: DoItVector has inverse order against GetPhysIntVector
    //       and SelectedAtRestDoItVector.
    //
    if(selectedAtRestDoItVector[fpProcessInfo->MAXofAtRestLoops - np - 1] != InActivated)
    {
      fpCurrentProcess =
          (G4VITProcess*) (*fpProcessInfo->fpAtRestDoItVector)[(G4int)np];

//      G4cout << " Invoke : "
//             << fpCurrentProcess->GetProcessName()
//             << G4endl;

      //  if(fpVerbose)
      //  {
      //    fpVerbose->AtRestDoItOneByOne();
      //  }

      fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess
          ->GetProcessID()));
      fpParticleChange = fpCurrentProcess->AtRestDoIt(*fpTrack, *fpStep);
      fpCurrentProcess->ResetProcessState();

      // Set the current process as a process which defined this Step length
      fpStep->GetPostStepPoint()->SetProcessDefinedStep(fpCurrentProcess);

      // Update Step
      fpParticleChange->UpdateStepForAtRest(fpStep);

      // Now Store the secondaries from ParticleChange to SecondaryList
      DealWithSecondaries(fN2ndariesAtRestDoIt);

      // Set the track status according to what the process defined
      // if kinetic energy >0, otherwise set  fStopButAlive
      fpTrack->SetTrackStatus(fpParticleChange->GetTrackStatus());

      // clear ParticleChange
      fpParticleChange->Clear();

    } //if(fSelectedAtRestDoItVector[np] != InActivated){
  } //for(std::size_t np=0; np < MAXofAtRestLoops; ++np){
  fpStep->UpdateTrack();

  // Modification par rapport au transport standard :
  // fStopAndKill doit etre propose par le modele
  // sinon d autres processus AtRest seront appeles
  // au pas suivant
  // fpTrack->SetTrackStatus(fStopAndKill);
}

//______________________________________________________________________________

// ************************************************************************
//  AlongStepDoIt
// ************************************************************************

void G4ITStepProcessor::InvokeAlongStepDoItProcs()
{

#ifdef DEBUG_MEM
  MemStat mem_first, mem_second, mem_diff;
#endif

#ifdef DEBUG_MEM
  mem_first = MemoryUsage();
#endif

  // If the current Step is defined by a 'ExclusivelyForced'
  // PostStepDoIt, then don't invoke any AlongStepDoIt
  if(fpState->fStepStatus == fExclusivelyForcedProc)
  {
    return;               // Take note 'return' at here !!!
  }

  // Invoke the all active continuous processes
  for(std::size_t ci = 0; ci < fpProcessInfo->MAXofAlongStepLoops; ++ci)
  {
    fpCurrentProcess =
        (G4VITProcess*) (*fpProcessInfo->fpAlongStepDoItVector)[(G4int)ci];
    if(fpCurrentProcess == 0) continue;
    // NULL means the process is inactivated by a user on fly.

    fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess
        ->GetProcessID()));
    fpParticleChange = fpCurrentProcess->AlongStepDoIt(*fpTrack, *fpStep);

#ifdef DEBUG_MEM
    MemStat mem_intermediaire = MemoryUsage();
    mem_diff = mem_intermediaire-mem_first;
    G4cout << "\t\t\t >> || MEM || After calling AlongStepDoIt for " << fpCurrentProcess->GetProcessName() << " and track "<< fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

//        fpCurrentProcess->SetProcessState(0);
    fpCurrentProcess->ResetProcessState();
    // Update the PostStepPoint of Step according to ParticleChange

    fpParticleChange->UpdateStepForAlongStep(fpStep);

#ifdef G4VERBOSE
    // !!!!! Verbose
    if(fpVerbose) fpVerbose->AlongStepDoItOneByOne();
#endif

    // Now Store the secondaries from ParticleChange to SecondaryList
    DealWithSecondaries(fN2ndariesAlongStepDoIt);

    // Set the track status according to what the process defined
    // if kinetic energy >0, otherwise set  fStopButAlive
    fpTrack->SetTrackStatus(fpParticleChange->GetTrackStatus());

    // clear ParticleChange
    fpParticleChange->Clear();
  }

#ifdef DEBUG_MEM
  MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire-mem_first;
  G4cout << "\t\t\t >> || MEM || After looping on processes with " << fpTrack->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

  fpStep->UpdateTrack();

  G4TrackStatus fNewStatus = fpTrack->GetTrackStatus();

  if(fNewStatus == fAlive && fpTrack->GetKineticEnergy() <= DBL_MIN)
  {
    //        G4cout << "G4ITStepProcessor::InvokeAlongStepDoItProcs : Track will be killed" << G4endl;
    if(fpProcessInfo->MAXofAtRestLoops>0) fNewStatus = fStopButAlive;
    else fNewStatus = fStopAndKill;
    fpTrack->SetTrackStatus( fNewStatus );
  }

}

//______________________________________________________________________________

// ************************************************************************
//  PostStepDoIt
// ************************************************************************

void G4ITStepProcessor::InvokePostStepDoItProcs()
{
  std::size_t _MAXofPostStepLoops = fpProcessInfo->MAXofPostStepLoops;
  G4SelectedPostStepDoItVector& selectedPostStepDoItVector = fpState
      ->fSelectedPostStepDoItVector;
  G4StepStatus& stepStatus = fpState->fStepStatus;

  // Invoke the specified discrete processes
  for(std::size_t np = 0; np < _MAXofPostStepLoops; ++np)
  {
    //
    // Note: DoItVector has inverse order against GetPhysIntVector
    //       and SelectedPostStepDoItVector.
    //
    G4int Cond = selectedPostStepDoItVector[_MAXofPostStepLoops
                                            - np - 1];
    if(Cond != InActivated)
    {
      if(((Cond == NotForced) && (stepStatus == fPostStepDoItProc)) ||
          ((Cond == Forced) && (stepStatus != fExclusivelyForcedProc))
         ||
         // ((Cond == Conditionally) && (stepStatus == fAlongStepDoItProc)) ||
         ((Cond == ExclusivelyForced) && (stepStatus == fExclusivelyForcedProc))
         || ((Cond == StronglyForced)))
      {

        InvokePSDIP(np);
      }
    } //if(*fSelectedPostStepDoItVector(np)........

    // Exit from PostStepLoop if the track has been killed,
    // but extra treatment for processes with Strongly Forced flag
    if(fpTrack->GetTrackStatus() == fStopAndKill)
    {
      for(std::size_t np1 = np + 1; np1 < _MAXofPostStepLoops; ++np1)
      {
        G4int Cond2 = selectedPostStepDoItVector[_MAXofPostStepLoops
                                                 - np1 - 1];
        if(Cond2 == StronglyForced)
        {
          InvokePSDIP(np1);
        }
      }
      break;
    }
  } //for(std::size_t np=0; np < MAXofPostStepLoops; ++np){
}

//______________________________________________________________________________

void G4ITStepProcessor::InvokePSDIP(std::size_t np)
{
  fpCurrentProcess = (G4VITProcess*) (*fpProcessInfo->fpPostStepDoItVector)[(G4int)np];

  fpCurrentProcess->SetProcessState(fpTrackingInfo->GetProcessState(fpCurrentProcess
      ->GetProcessID()));
  fpParticleChange = fpCurrentProcess->PostStepDoIt(*fpTrack, *fpStep);
//    fpCurrentProcess->SetProcessState(0);
  fpCurrentProcess->ResetProcessState();

  // Update PostStepPoint of Step according to ParticleChange
  fpParticleChange->UpdateStepForPostStep(fpStep);

#ifdef G4VERBOSE
  // !!!!! Verbose
  if(fpVerbose) fpVerbose->PostStepDoItOneByOne();
#endif

  // Update G4Track according to ParticleChange after each PostStepDoIt
  fpStep->UpdateTrack();

  // Update safety after each invocation of PostStepDoIts
  fpStep->GetPostStepPoint()->SetSafety(CalculateSafety());

  // Now Store the secondaries from ParticleChange to SecondaryList
  DealWithSecondaries(fN2ndariesPostStepDoIt);

  // Set the track status according to what the process defined
  fpTrack->SetTrackStatus(fpParticleChange->GetTrackStatus());

  // clear ParticleChange
  fpParticleChange->Clear();
}

//______________________________________________________________________________

// ************************************************************************
//  Transport on a given time
// ************************************************************************

void G4ITStepProcessor::FindTransportationStep()
{
  double physicalStep(0.);

  fpTransportation = fpProcessInfo->fpTransportation;
  // dynamic_cast<G4ITTransportation*>((fpProcessInfo->fpAlongStepGetPhysIntVector)[MAXofAlongStepLoops-1]);

  if(!fpTrack)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No G4ITStepProcessor::fpTrack found";
    G4Exception("G4ITStepProcessor::FindTransportationStep",
                "ITStepProcessor0013",
                FatalErrorInArgument,
                exceptionDescription);
    return;

  }
  if(!fpITrack)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No G4ITStepProcessor::fITrack";
    G4Exception("G4ITStepProcessor::FindTransportationStep",
                "ITStepProcessor0014",
                FatalErrorInArgument,
                exceptionDescription);
    return;
  }
  if(!(fpITrack->GetTrack()))
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "No G4ITStepProcessor::fITrack->GetTrack()";
    G4Exception("G4ITStepProcessor::FindTransportationStep",
                "ITStepProcessor0015",
                FatalErrorInArgument,
                exceptionDescription);
    return;
  }

  if(fpTransportation)
  {
    fpTransportation->SetProcessState(fpTrackingInfo->GetProcessState(fpTransportation
        ->GetProcessID()));
    fpTransportation->ComputeStep(*fpTrack, *fpStep, fTimeStep, physicalStep);

//        fpTransportation->SetProcessState(0);
    fpTransportation->ResetProcessState();
  }

  if(physicalStep >= DBL_MAX)
  {
//    G4cout << "---- 2) Setting stop and kill for " << GetIT(fpTrack)->GetName() << G4endl;
    fpTrack -> SetTrackStatus(fStopAndKill);
    return;
  }

  fpState->fPhysicalStep = physicalStep;
}

//______________________________________________________________________________

void G4ITStepProcessor::InvokeTransportationProc()
{
  std::size_t _MAXofPostStepLoops = fpProcessInfo->MAXofPostStepLoops;
  G4SelectedPostStepDoItVector& selectedPostStepDoItVector = fpState
      ->fSelectedPostStepDoItVector;
  G4StepStatus& stepStatus = fpState->fStepStatus;

  // Invoke the specified discrete processes
  for(std::size_t np = 0; np < _MAXofPostStepLoops; ++np)
  {
    //
    // Note: DoItVector has inverse order against GetPhysIntVector
    //       and SelectedPostStepDoItVector.
    //
    G4int Cond = selectedPostStepDoItVector[_MAXofPostStepLoops - np - 1];
    if(Cond != InActivated)
    {
      if(((Cond == Forced) && (stepStatus != fExclusivelyForcedProc)) ||
      // ((Cond == Conditionally) && (stepStatus == fAlongStepDoItProc)) ||
      ((Cond == ExclusivelyForced) && (stepStatus == fExclusivelyForcedProc))
         || ((Cond == StronglyForced)))
      {

        InvokePSDIP(np);
      }
    } //if(Cond != InActivated)

    // Exit from PostStepLoop if the track has been killed,
    // but extra treatment for processes with Strongly Forced flag
    if(fpTrack->GetTrackStatus() == fStopAndKill)
    {
      for(std::size_t np1 = np + 1; np1 < _MAXofPostStepLoops; ++np1)
      {
        G4int Cond2 = selectedPostStepDoItVector[_MAXofPostStepLoops - np1 - 1];
        if(Cond2 == StronglyForced)
        {
          InvokePSDIP(np1);
        }
      }
      break;
    }
  }
}

//______________________________________________________________________________

// ************************************************************************
//  Apply cuts
// ************************************************************************

void G4ITStepProcessor::ApplyProductionCut(G4Track* aSecondary)
{
  G4bool tBelowCutEnergyAndSafety = false;
  G4int tPtclIdx = G4ProductionCuts::GetIndex(aSecondary->GetDefinition());
  if(tPtclIdx < 0)
  {
    return;
  }
  G4ProductionCutsTable* tCutsTbl =
      G4ProductionCutsTable::GetProductionCutsTable();
  G4int tCoupleIdx = tCutsTbl->GetCoupleIndex(fpPreStepPoint
      ->GetMaterialCutsCouple());
  G4double tProdThreshold =
      (*(tCutsTbl->GetEnergyCutsVector(tPtclIdx)))[tCoupleIdx];
  if(aSecondary->GetKineticEnergy() < tProdThreshold)
  {
    tBelowCutEnergyAndSafety = true;
    if(std::abs(aSecondary->GetDynamicParticle()->GetCharge()) > DBL_MIN)
    {
      G4double currentRange
      = G4LossTableManager::Instance()->GetRange(aSecondary->GetDefinition(),
          aSecondary->GetKineticEnergy(),
          fpPreStepPoint->GetMaterialCutsCouple());
      tBelowCutEnergyAndSafety = (currentRange < CalculateSafety() );
    }
  }

  if(tBelowCutEnergyAndSafety)
  {
    if(!(aSecondary->IsGoodForTracking()))
    {
      // Add kinetic energy to the total energy deposit
      fpStep->AddTotalEnergyDeposit(aSecondary->GetKineticEnergy());
      aSecondary->SetKineticEnergy(0.0);
    }
  }
}
