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
// $Id: G4Scheduler.cc 60494 2012-07-12 14:49:30Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include <G4AllITFinder.hh>
#include <G4Scheduler.hh>
#include <G4SchedulerMessenger.hh>

#include "G4SystemOfUnits.hh"
#include "G4ITModelProcessor.hh"
#include "G4ITStepProcessor.hh"
#include "G4IT.hh"
#include "G4ITReactionChange.hh"
#include "G4ITModelHandler.hh"
#include "G4VITStepModel.hh"
#include "G4UserTimeStepAction.hh"
#include "G4ITTrackingManager.hh"
#include "G4ITTrackingInteractivity.hh"
#include "G4TrackingInformation.hh"
#include "G4UnitsTable.hh"
#include "G4ITStepStatus.hh"
#include "G4ITGun.hh"
#include "G4StateManager.hh"
#include "G4Timer.hh"
#include <sstream>

//#define DEBUG_MEM 1
#define DEBUG_MEM_STEPPING
#define DEBUG_MEM_DETAILED_STEPPING
//#define DEBUG_MEM_DOIT

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
using G4MemStat::MemStat;
#endif

//COLOR FOR DEBUGING
//#define USE_COLOR 1

#ifdef USE_COLOR
#define RED  "\033[0;31m"
#define LIGHT_RED  "\33[1;31m"
#define GREEN "\033[32;40m"
#define GREEN_ON_BLUE "\033[1;32;44m"
#define RESET_COLOR "\033[0m"
#else
#define RED  ""
#define LIGHT_RED  ""
#define GREEN ""
#define GREEN_ON_BLUE ""
#define RESET_COLOR ""
#endif

using namespace std;

G4ThreadLocal G4Scheduler* G4Scheduler::fgScheduler(0);

class IosFlagSaver
{
public:
  explicit IosFlagSaver(std::ostream& _ios) :
      ios(_ios), f(_ios.flags())
  {
  }
  ~IosFlagSaver()
  {
    ios.flags(f);
  }

  // IosFlagSaver(const IosFlagSaver &rhs) = delete;
  // IosFlagSaver& operator= (const IosFlagSaver& rhs) = delete;

private:
  std::ostream& ios;
  std::ios::fmtflags f;
};

template<typename T>
  inline bool IsInf(T value)
  {
    return std::numeric_limits<T>::has_infinity
        && value == std::numeric_limits<T>::infinity();
  }
//_________________________________________________________________________

G4Scheduler* G4Scheduler::Instance()
{
  if (fgScheduler == 0) fgScheduler = new G4Scheduler();
  return fgScheduler;
}
//_________________________________________________________________________

G4bool G4Scheduler::Notify(G4ApplicationState requestedState)
{
  if (requestedState == G4State_Quit)
  {
    if (fVerbose >= 4)
    {
      G4cout << "G4Scheduler received G4State_Quit" << G4endl;
    }
    Clear();
    //DeleteInstance();
  }
  return true;
}
//_________________________________________________________________________

void G4Scheduler::DeleteInstance()
{
  if (fgScheduler)
  {
    delete fgScheduler;
  }
}
//_________________________________________________________________________

G4Scheduler::G4Scheduler() :
    G4VScheduler(), G4VStateDependent(),
    fTrackContainer((G4ITTrackHolder&)*G4ITTrackHolder::Instance())
{
  Create();
}

void G4Scheduler::Create()
{
  fUseDefaultTimeSteps = true;
  fUserUpperTimeLimit = -1;
  fpGun = 0;
  fContinue = true;
//  fpMainList = 0;
//  fpWaitingList = 0;
  fpTrackingInteractivity = 0;

  fITStepStatus = eUndefined;

  fpUserTimeSteps = 0;

  fTimeStep = DBL_MAX;
  fTSTimeStep = DBL_MAX;
  fILTimeStep = DBL_MAX;
  fPreviousTimeStep = DBL_MAX;

  fZeroTimeCount = 0;
  fMaxNZeroTimeStepsAllowed = 10;

  fStartTime = 0;
  fComputeTimeStep = false;
  fComputeReaction = false;
  fTimeTolerance = 1 * picosecond;
  fEndTime = 1 * microsecond;
  fGlobalTime = -1;
  fInteractionStep = true;
  fUsePreDefinedTimeSteps = false;

  fDefaultMinTimeStep = 1 * picosecond;

  fpStepProcessor = 0;
  //    fpMasterStepProcessor = 0 ;

  fpModelProcessor = 0;
  //    fpMasterModelProcessor= 0;

  fNbSteps = 0;
  fMaxSteps = -1;

  fRunning = false;
  fInitialized = false;

  fpUserTimeStepAction = 0;
  fpTrackingManager = 0;

  fNbTracks = -1;
  fpModelHandler = new G4ITModelHandler();
  fpTrackingManager = new G4ITTrackingManager();

  fVerbose = 0;
  fWhyDoYouStop = false;
  fDefinedMinTimeStep = -1.;
  fReachedUserTimeLimit = false;
  fTmpEndTime = -1.;
  fTmpGlobalTime = -1.;

  fSteppingMsg = new G4SchedulerMessenger(this);

  G4ITTypeManager::Instance()->ReserveRessource();
}

//_________________________________________________________________________

G4Scheduler::~G4Scheduler()
{

  if (fSteppingMsg) // is used as a flag to know whether the manager was cleared
  {
    Clear();
  }
  fgScheduler = 0;

//  if (fVerbose >= 1)
//  {
//    G4cout << "G4Scheduler is being deleted. Bye :) !" << G4endl;
//  }
}

void G4Scheduler::Clear()
{
//  if (fVerbose) G4cout << "*** G4Scheduler is being cleared ***" << G4endl;

  if(fSteppingMsg)
  {
    delete fSteppingMsg;
    fSteppingMsg = 0;
  }
  if(fpStepProcessor)
  {
    delete fpStepProcessor;
    fpStepProcessor = 0;
  }
  if(fpModelProcessor)
  {
    delete fpModelProcessor;
    fpModelProcessor = 0;
  }
  if(fpModelHandler)
  {
    delete fpModelHandler;
    fpModelHandler = 0;
  }
  //    if(fpMasterStepProcessor) delete fpMasterStepProcessor ;
  //    if(fpMasterModelProcessor) delete fpMasterModelProcessor ;
  G4ITTypeManager::Instance() -> ReleaseRessource();
  ClearList();
  if(fpTrackingManager)
  {
    delete fpTrackingManager;
    fpTrackingManager = 0;
  }

  /*
   * DEBUG
   * assert(G4StateManager::GetStateManager()->
   * DeregisterDependent(this) == true);
   */
}

//_________________________________________________________________________

void G4Scheduler::ClearList()
{
//  if (fNbTracks == 0) return;

  fTrackContainer.Clear();

  fNbTracks = -1;

  G4AllITFinder::DeleteInstance();
}

//_________________________________________________________________________
void G4Scheduler::RegisterModel(G4VITStepModel* model, double time)
{
  fpModelHandler->RegisterModel(model, time);
}

//_________________________________________________________________________

void G4Scheduler::Initialize()
{
  if (fpStepProcessor)
  {
    delete fpStepProcessor;
  }
  if (fpModelProcessor)
  {
    delete fpModelProcessor;
  }
  // if(fpMasterModelProcessor)
  // {
  //     delete fpMasterModelProcessor;
  // }

  //______________________________________________________________

  fpModelProcessor = new G4ITModelProcessor();
  fpModelProcessor->SetModelHandler(fpModelHandler);

  // fpMasterModelProcessor = new G4ITModelProcessor();
  // fpMasterModelProcessor->SetModelHandler(fpModelHandler);
  // fpModelProcessor      = fpMasterModelProcessor;

  //______________________________________________________________

  fpStepProcessor = new G4ITStepProcessor();
  fpStepProcessor->SetTrackingManager(fpTrackingManager);

  fpTrackingManager->SetInteractivity(fpTrackingInteractivity);

  // fpMasterStepProcessor = new G4ITStepProcessor();
  // fpMasterStepProcessor->SetTrackingManager(fpTrackingManager);
  // fpStepProcessor = fpMasterStepProcessor ;

  if (fpModelHandler->GetTimeStepComputerFlag()) fComputeTimeStep = true;
  if (fpModelHandler->GetReactionProcessFlag()) fComputeReaction = true;
  //______________________________________________________________

  fInitialized = true;
}
//_________________________________________________________________________

void G4Scheduler::Reset()
{
  fStartTime = 0;
  fUserUpperTimeLimit = -1;
  fTimeStep = DBL_MAX;
  fTSTimeStep = DBL_MAX;
  fILTimeStep = DBL_MAX;
  fPreviousTimeStep = DBL_MAX;
  fGlobalTime = -1;
  fInteractionStep = true;
  fITStepStatus = eUndefined;
  fZeroTimeCount = 0;

  fNbSteps = 0;
  fContinue = true;
  // fReactingTracks.clear();
  fLeadingTracks.clear();
  fReactionSet.CleanAllReaction();
}
//_________________________________________________________________________

void G4Scheduler::Process()
{

#ifdef G4VERBOSE
  if (fVerbose)
  {
    G4cout << "*** G4Scheduler starts processing " << G4endl;
    if(fVerbose > 2)
    G4cout << "___________________________________________"
    "___________________________" << G4endl;
  }
#endif

  if (fInitialized == false) Initialize();

  // fpTrackingManager->Initialize();
  fpModelProcessor->Initialize();
  fpStepProcessor->Initialize();

  // TODO
  // fpMasterModelProcessor->Initialize();
  // fpMasterStepProcessor->Initialize();

  if (fpGun) fpGun->DefineTracks();

  if (fpTrackingInteractivity) fpTrackingInteractivity->Initialize();

  // ___________________
  fRunning = true;
  Reset();

  if (fpUserTimeStepAction) fpUserTimeStepAction->StartProcessing();

#ifdef G4VERBOSE
  G4bool trackFound = false;
#endif

  //===========================================================================
  // By default, before the G4Scheduler is launched, the tracks are pushed to
  // the delayed lists
  //===========================================================================

  if(fTrackContainer.DelayListsNOTEmpty())
  {
    fStartTime = fTrackContainer.GetNextTime();
#ifdef G4VERBOSE
    trackFound = true;
    G4Timer localtimer;
    if(fVerbose>1)
    {
      localtimer.Start();
    }
#endif
    SynchronizeTracks();
#ifdef G4VERBOSE
    if(fVerbose>1)
    {
      localtimer.Stop();
      G4cout << "G4Scheduler: process time= "<< localtimer << G4endl;
    }
#endif
  }

//  //---------------------------------
//  // TODO: This could be removed ?
//  if(fTrackContainer.MainListsNOTEmpty()
//      && (fMaxSteps == -1 ? true : fNbSteps < fMaxSteps)
//      && fGlobalTime < fEndTime
//      && fContinue == true)
//{
//#ifdef G4VERBOSE
//    trackFound = true;
//#endif
//    DoProcess();
//}
//  //---------------------------------

#ifdef G4VERBOSE
  if(fVerbose)
  {
    if(trackFound)
    {
      G4cout << "*** G4Scheduler ends at time : "
      << G4BestUnit(fGlobalTime,"Time") << G4endl;
      G4cout << "___________________________________" << G4endl;
    }
    else
    {
      G4cout << "*** G4Scheduler did not start because no "
      "track was found to be processed"<< G4endl;
      G4cout << "___________________________________" << G4endl;
    }
  }
#endif

  fRunning = false;

  if (fpUserTimeStepAction) fpUserTimeStepAction->EndProcessing();

  // ___________________
  EndTracking();
  ClearList();

  Reset();

  if (fpTrackingInteractivity) fpTrackingInteractivity->Finalize();
}
//_________________________________________________________________________

void G4Scheduler::SynchronizeTracks()
{
//  if(fTrackContainer.WaitingListsNOTEmpty())
//   {
//    G4ExceptionDescription exceptionDescription;
//    exceptionDescription
//        << "There is a waiting track list (fpWaitingList != 0).";
//    exceptionDescription
//        << " When G4Scheduler::SynchronizeTracks() is called, ";
//    exceptionDescription
//        << "no more tracks should remain in the fpWaitingList.";
//    G4Exception("G4Scheduler::SynchronizeTracks", "ITScheduler002",
//                FatalErrorInArgument, exceptionDescription);
//  }

  // Backup main list and time feature
  // Reminder : the main list here, should
  // have the biggest global time
  //  fTrackContainer.MoveMainToWaitingList();
  // TODO: not yet supported

  fTmpGlobalTime = fGlobalTime;
  fTmpEndTime = fEndTime;

  fGlobalTime = fTrackContainer.GetNextTime();
  G4double tmpGlobalTime = fGlobalTime;
  while(fTrackContainer.MergeNextTimeToMainList(tmpGlobalTime))
  {
//    assert(tmpGlobalTime == fGlobalTime);
    fEndTime =  min(fTrackContainer.GetNextTime(), fTmpEndTime);
    DoProcess();
  }
}
//_________________________________________________________________________

void G4Scheduler::DoProcess()
// We split it from the Process() method to avoid repeating code in SynchronizeTracks
{
  if (fpUserTimeStepAction) fpUserTimeStepAction->NewStage();

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_STEPPING)
  MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_STEPPING)
  mem_first = MemoryUsage();
#endif

  while (fGlobalTime < fEndTime
         && fTrackContainer.MainListsNOTEmpty()
         && (fMaxSteps == -1 ? true : fNbSteps < fMaxSteps)
         && fContinue == true)
  {
//    G4cout << "Mainlist size : " << fTrackContainer.GetMainList()->size()
//        << G4endl;

    Stepping();

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_STEPPING)
    mem_second = MemoryUsage();
    mem_diff = mem_second-mem_first;
    G4cout << "\t || MEM || After step " << fNbSteps << ", diff is : "
    << mem_diff << G4endl;
#endif
  }

#ifdef G4VERBOSE
  if (fWhyDoYouStop)
  {
    G4cout << "G4Scheduler has reached a stage: it might be"
           " a transition or the end"
           << G4endl;

    G4bool normalStop = false;

    if(fGlobalTime >= fEndTime)
    {
      G4cout << "== G4Scheduler: I stop because I reached the final time : "
      << G4BestUnit(fEndTime,"Time") << " =="<< G4endl;
      normalStop = true;
    }
    if(fTrackContainer.MainListsNOTEmpty()  == false) // is empty
    {
      G4cout << "G4Scheduler: I stop because the current main list of tracks"
      "is empty"
      << G4endl;
      normalStop = true;
    }
    if(fMaxSteps == -1 ? false : fNbSteps >= fMaxSteps)
    {
      G4cout << "G4Scheduler: I stop because I reached the maximum allowed "
      "number of steps=" << fMaxSteps
      << G4endl;
      normalStop = true;
    }
    if(fContinue && normalStop == false)
    {
      G4cout << "G4Scheduler: It might be that I stop because "
      "I have been told so. You may check "
      "member fContinue and usage of the method G4Scheduler::Stop()."
      << G4endl;
    }
  }
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "\t || MEM || After stepping, diff is : " << mem_diff << G4endl;
#endif

#ifdef G4VERBOSE
  if (fVerbose > 2) G4cout
      << "*** G4Scheduler has finished processing a track list at time : "
      << G4BestUnit(fGlobalTime, "Time") << G4endl;
#endif
}
//_________________________________________________________________________

void G4Scheduler::Stepping()
{
  IosFlagSaver iosfs(G4cout);
  fTimeStep = DBL_MAX;

  fTSTimeStep = DBL_MAX;
  fILTimeStep = DBL_MAX;

  fInteractionStep = false;
  fReachedUserTimeLimit = false;

  fITStepStatus = eUndefined;

  // Start of step
#ifdef G4VERBOSE
  if (fVerbose > 2)
  {
#ifdef USE_COLOR
    G4cout << LIGHT_RED;
#endif
    G4cout << "*** Start Of Step N°" << fNbSteps + 1 << " ***" << G4endl;
    G4cout << "Current Global time : " << G4BestUnit(fGlobalTime, "Time")
           <<G4endl;
#ifdef USE_COLOR
    G4cout << RESET_COLOR;
#endif
  }
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  fDefinedMinTimeStep = GetLimitingTimeStep();

  if (fUsePreDefinedTimeSteps)
  {
    if (fpUserTimeSteps == 0) // Extra checking, is it necessary ?
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription
        << "You are asking to use user defined steps but you did not give any.";
      G4Exception("G4Scheduler::FindUserPreDefinedTimeStep",
                  "ITScheduler004", FatalErrorInArgument,
                  exceptionDescription);
      return; // makes coverity happy
    }
#ifdef G4VERBOSE
    if (fVerbose > 2)
    {
#ifdef USE_COLOR
      G4cout << LIGHT_RED;
#endif
      G4cout << "*** At time : " << G4BestUnit(fGlobalTime, "Time")
             << " the chosen user time step is : "
             << G4BestUnit(fDefinedMinTimeStep, "Time") << " ***" << G4endl;
#ifdef USE_COLOR
      G4cout << RESET_COLOR;
#endif
    }
#endif
  }

  if (fComputeTimeStep)
  {
    CalculateMinTimeStep(); // => at least N (N = nb of tracks) loops
  }
  else if(fUseDefaultTimeSteps)
  {
    fTSTimeStep = fDefinedMinTimeStep;
  }

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "|| MEM || After computing TS, diff is : " << mem_diff << G4endl;
#endif

#ifdef G4VERBOSE
  if (fVerbose > 2)
  {
#ifdef USE_COLOR
    G4cout << LIGHT_RED;
#endif
    G4cout << "*** Time stepper returned : " << G4BestUnit(fTSTimeStep, "Time")
           << " ***" << G4endl;
#ifdef USE_COLOR
    G4cout << RESET_COLOR;
#endif
  }
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  // Call IL even if fTSTimeStep == 0
  // if fILTimeStep == 0 give the priority to DoIt processes
  ComputeInteractionLength(); // => at least N loops
  // All process returns the physical step of interaction
  // using the approximation : E = cste along the step
  // Only the transportation calculates the corresponding
  // time step

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "|| MEM || After IL, diff is : " << mem_diff << G4endl;
#endif

#ifdef G4VERBOSE
  if (fVerbose > 2)
  {
#ifdef USE_COLOR
    G4cout << LIGHT_RED;
#endif
    G4cout << "*** The minimum time returned by the processes is : "
           << G4BestUnit(fILTimeStep, "Time") << " ***" << G4endl;
#ifdef USE_COLOR
    G4cout << RESET_COLOR;
#endif
  }
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  if (fILTimeStep <= fTSTimeStep)
    // Give the priority to the IL
  {
    fInteractionStep = true;
//    fReactingTracks.clear(); // Give the priority to the IL
    fReactionSet.CleanAllReaction();
    fTimeStep = fILTimeStep;
    fITStepStatus = eInteractionWithMedium;
  }
  else
  {
    fInteractionStep = false;
    ResetLeadingTracks();
    fTimeStep = fTSTimeStep;
    fITStepStatus = eCollisionBetweenTracks;
  }

  if (fGlobalTime + fTimeStep > fEndTime)
  // This check is done at very time step
  {
    fTimeStep = fEndTime - fGlobalTime;
    fITStepStatus = eInteractionWithMedium; // ie: transportation
    fInteractionStep = true;
//    fReactingTracks.clear();
    fReactionSet.CleanAllReaction();
    ResetLeadingTracks();
  }

  if (fTimeStep == 0) // < fTimeTolerance)
  {
    fZeroTimeCount++;
    if (fZeroTimeCount >= fMaxNZeroTimeStepsAllowed)
    {
      G4ExceptionDescription exceptionDescription;

      exceptionDescription << "Too many zero time steps were detected. ";
      exceptionDescription << "The simulation is probably stuck. ";
      exceptionDescription
          << "The maximum number of zero time steps is currently : "
          << fMaxNZeroTimeStepsAllowed;
      exceptionDescription << ".";

      G4Exception("G4Scheduler::Stepping", "ITSchedulerNullTimeSteps",
                  FatalErrorInArgument, exceptionDescription);
    }
  }
  else
  {
    fZeroTimeCount = 0;
  }

  fReachedUserTimeLimit =
      ((fTimeStep <= fDefinedMinTimeStep) || ((fTimeStep > fDefinedMinTimeStep)
          && fabs(fTimeStep - fDefinedMinTimeStep) < fTimeTolerance)) ?
          true : false;

  if (fpUserTimeStepAction) fpUserTimeStepAction->UserPreTimeStepAction();
  // TODO: pre/post

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "|| MEM || After LeadingTracks and UserPreTimeStepAction: "
         << mem_diff << G4endl;
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  // if fTSTimeStep > 0 => needs to call the transportation process
  // if fILTimeStep < fTSTimeStep => call DoIt processes
  // if fILTimeStep == fTSTimeStep => give the priority to the DoIt processes
  if (fTSTimeStep > 0 || fILTimeStep <= fTSTimeStep)
  {
    //        G4cout << "Will call DoIT" << G4endl;
    DoIt();

    fTrackContainer.MergeSecondariesWithMainList();
    KillTracks();
  }
  // else
  // {
  //   G4cout << "fTSTimeStep : " << fTSTimeStep << G4endl;
  //   G4cout << "fILTimeStep : " << fILTimeStep << G4endl;
  // }

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "|| MEM || After DoIT, diff is : " << mem_diff << G4endl;
#endif

  fGlobalTime += fTimeStep;

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  ComputeTrackReaction();

  fTrackContainer.MergeSecondariesWithMainList();

  fTrackContainer.KillTracks();

  fNbSteps++;

  if (fpUserTimeStepAction)
  {
    fpUserTimeStepAction->UserPostTimeStepAction();
  }

  fPreviousTimeStep = fTimeStep;

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "|| MEM || After computing reactions + UserPostTimeStepAction, "
  "diff is : " << mem_diff << G4endl;
#endif

  // End of step
#ifdef G4VERBOSE
    if (fVerbose >= 2)
    {
#ifdef USE_COLOR
      G4cout << LIGHT_RED;
#endif

      G4String interactionType;
      GetCollisionType(interactionType);

      std::stringstream finalOutput;

      finalOutput << "*** End of step N°" << fNbSteps
             << "\t T_i= " << G4BestUnit(fGlobalTime-fTimeStep, "Time")
             << "\t dt= " << G4BestUnit(fTimeStep, "Time")
             << "\t T_f= " << G4BestUnit(fGlobalTime, "Time")
             << "\t " << interactionType
             << G4endl;

      if(fVerbose>2)
      {
        if(fReachedUserTimeLimit)
        {
          finalOutput << "It has also reached the user time limit" << G4endl;
        }
        finalOutput << "_______________________________________________________________"
            "_______"<< G4endl;
      }

      G4cout << finalOutput.str();

#ifdef USE_COLOR
      G4cout << RESET_COLOR;
#endif
    }
#endif

}
//_________________________________________________________________________

double G4Scheduler::GetLimitingTimeStep() const
{
  if (fpUserTimeSteps == 0) return fDefaultMinTimeStep;
  if (fabs(fGlobalTime - fUserUpperTimeLimit) < fTimeTolerance)
    return fDefinedMinTimeStep;

  map<double, double>::const_iterator it_fpUserTimeSteps_i = fpUserTimeSteps
      ->upper_bound(fGlobalTime);
  map<double, double>::const_iterator it_fpUserTimeSteps_low = fpUserTimeSteps
      ->lower_bound(fGlobalTime);

  // DEBUG
  // G4cout << "fGlobalTime : " << G4BestUnit(fGlobalTime,"Time")
  //            << G4endl;
  // G4cout << "fpUserTimeSteps_i : "
  //        <<"<"<<G4BestUnit(it_fpUserTimeSteps->first,"Time")
  //        <<", "<< G4BestUnit(it_fpUserTimeSteps->second,"Time")<<">"
  //        << "\t fpUserTimeSteps_low : "
  //        <<"<"<<G4BestUnit(fpUserTimeSteps_low->first,"Time")<<", "*
  //        << G4BestUnit(fpUserTimeSteps_low->second,"Time")<<">"
  //        << G4endl;

  if (it_fpUserTimeSteps_i == fpUserTimeSteps->end())
  {
    it_fpUserTimeSteps_i--;
    fUserUpperTimeLimit = fEndTime;
  }
  else if (fabs(fGlobalTime - it_fpUserTimeSteps_low->first) < fTimeTolerance)
  {
    // Case : fGlobalTime = X picosecond
    // and fpUserTimeSteps_low->first = X picosecond
    // but the precision is not good enough
    it_fpUserTimeSteps_i = it_fpUserTimeSteps_low;
    map<double, double>::const_iterator tmp_it = it_fpUserTimeSteps_low;
    tmp_it++;
    if (tmp_it == fpUserTimeSteps->end())
    {
      fUserUpperTimeLimit = fEndTime;
    }
    else
    {
      fUserUpperTimeLimit = tmp_it->second;
    }
  }
  else if (it_fpUserTimeSteps_i == it_fpUserTimeSteps_low)
  {
    // "Normal" cases
    fUserUpperTimeLimit = it_fpUserTimeSteps_i->second;
    it_fpUserTimeSteps_i--;
  }
  else
  {
    fUserUpperTimeLimit = it_fpUserTimeSteps_i->second;
    it_fpUserTimeSteps_i = it_fpUserTimeSteps_low;
  }

  return it_fpUserTimeSteps_i->second;
}

//_________________________________________________________________________

void G4Scheduler::FindUserPreDefinedTimeStep()
{

  if (fpUserTimeSteps == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "You are asking to use user defined steps but you did not give any.";
    G4Exception("G4Scheduler::FindUserPreDefinedTimeStep",
                "ITScheduler004", FatalErrorInArgument, exceptionDescription);
    return; // makes coverity happy
  }
  map<double, double>::iterator fpUserTimeSteps_i =
      fpUserTimeSteps->upper_bound(fGlobalTime);
  map<double, double>::iterator fpUserTimeSteps_low = fpUserTimeSteps
      ->lower_bound(fGlobalTime);

  // DEBUG
  // G4cout << "fGlobalTime : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
  // G4cout << "fpUserTimeSteps_i : "
  // <<"<"<<G4BestUnit(fpUserTimeSteps_i->first,"Time")<<", "
  // << G4BestUnit(fpUserTimeSteps_i->second,"Time")<<">"
  // << "\t fpUserTimeSteps_low : "
  // <<"<"<<G4BestUnit(fpUserTimeSteps_low->first,"Time")<<", "
  // << G4BestUnit(fpUserTimeSteps_low->second,"Time")<<">"
  // << G4endl;

  if (fpUserTimeSteps_i == fpUserTimeSteps->end())
  {
    fpUserTimeSteps_i--;
  }
  else if (fabs(fGlobalTime - fpUserTimeSteps_low->first) < fTimeTolerance)
  {
    // Case : fGlobalTime = X picosecond
    // and fpUserTimeSteps_low->first = X picosecond
    // but the precision is not good enough
    fpUserTimeSteps_i = fpUserTimeSteps_low;
  }
  else if (fpUserTimeSteps_i == fpUserTimeSteps_low)
  {
    // "Normal" cases
    fpUserTimeSteps_i--;
  }
  else
  {
    fpUserTimeSteps_i = fpUserTimeSteps_low;
  }

  fDefinedMinTimeStep = fpUserTimeSteps_i->second;
}
//_________________________________________________________________________

void G4Scheduler::CalculateMinTimeStep()
{

  if (fpModelProcessor == 0)
  //        if(fpMasterModelProcessor == 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "There is no G4ITModelProcessor to hande IT reaction. ";
    exceptionDescription
        << "You probably did not initialize the G4Scheduler. ";
    exceptionDescription
        << "Just do G4Scheduler::Instance()->Initialize(); ";
    exceptionDescription << " but only after initializing the run manager.";
    G4Exception("G4Scheduler::CalculateMinStep", "ITScheduler005",
                FatalErrorInArgument, exceptionDescription);
    return; // makes coverity happy
  }

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_first = MemoryUsage();
#endif

  fpModelProcessor->InitializeStepper(fGlobalTime, fDefinedMinTimeStep);
  // TODO
  // fpMasterModelProcessor -> InitializeStepper(fGlobalTime, fDefinedMinTimeStep) ;

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "\t || MEM || G4Scheduler::CalculateMinTimeStep || After "
      "computing fpModelProcessor -> InitializeStepper, diff is : "
      << mem_diff
      << G4endl;
#endif

//  G4TrackList::iterator fpMainList_i = fpMainList->begin();

  G4TrackManyList* mainList = fTrackContainer.GetMainList();
  G4TrackManyList::iterator it = mainList->begin();
  G4TrackManyList::iterator end = mainList->end();

  for (; it != end; it++)
  {
    G4Track * track = *it;

    if (track == 0)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No track found.";
      G4Exception("G4Scheduler::CalculateMinStep", "ITScheduler006",
                  FatalErrorInArgument, exceptionDescription);
      return; // makes coverity happy
    }

#ifdef DEBUG
    G4cout << "*_* " << GetIT(track)->GetName()
        << " ID: " << track->GetTrackID()
        <<  " at time : " << track->GetGlobalTime()
        << G4endl;
#endif

    G4TrackStatus trackStatus = track->GetTrackStatus();
    if (trackStatus == fStopAndKill || trackStatus == fStopButAlive)
    {
      continue;
    }

    // This Extract... was thought for MT mode at the track level
    //ExtractTimeStepperData(fpModelProcessor) ;

    fpModelProcessor->CalculateTimeStep(track, fDefinedMinTimeStep);

    // if MT mode at track level, this command should be displaced
    ExtractTimeStepperData(fpModelProcessor);
  }

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DETAILED_STEPPING)
  mem_second = MemoryUsage();
  mem_diff = mem_second-mem_first;
  G4cout << "\t || MEM || G4Scheduler::CalculateMinTimeStep || "
      "After looping on tracks, diff is : " << mem_diff << G4endl;
#endif
}
//_________________________________________________________________________

void G4Scheduler::ExtractTimeStepperData(G4ITModelProcessor* MP)
{
  G4Track* track = (G4Track*) MP->GetTrack();
  if (track == 0)
  {
    MP->CleanProcessor();
    return;
  }

  const std::vector<std::vector<G4VITStepModel*> >* model =
      MP->GetCurrentModel();

  for (unsigned i = 0; i < model->size(); i++)
  {
    for (unsigned j = 0; j < (*model)[i].size(); j++)
    {
      G4VITStepModel* mod = (*model)[i][j];

      if (mod == 0)
      {
        continue;
      }

      G4VITTimeStepComputer* stepper(mod->GetTimeStepper());

      G4double sampledMinTimeStep(stepper->GetSampledMinTimeStep());
      G4TrackVectorHandle reactants(stepper->GetReactants());

      if (sampledMinTimeStep < fTSTimeStep)
      {
        /*
         // DEBUG SPECIAL CASE
         if(!reactants)
         {
         G4ExceptionDescription exceptionDescription ;
         exceptionDescription << "No reactants were found by the time stepper.";
         G4Exception("G4Scheduler::ExtractTimeStepperData","ITScheduler007",
         FatalErrorInArgument,exceptionDescription);
         continue ;
         }
         */

        fTSTimeStep = sampledMinTimeStep;
        //fReactingTracks.clear();
        fReactionSet.CleanAllReaction();
        if (bool(reactants))
        {
          // G4cout << "*** (1) G4Scheduler::ExtractTimeStepperData insert
          //          reactants for " << GetIT(track)->GetName() << G4endl;
          // G4cout << "bool(reactants) = " << bool(reactants) << G4endl;
          // G4cout << reactants->size() << G4endl;
          // G4cout << GetIT(reactants->operator[](0))->GetName() << G4endl;

         //  fReactingTracks.insert(make_pair(track, reactants));
          fReactionSet.AddReactions(fTSTimeStep, track, reactants);
          stepper->ResetReactants();
        }
      }
      else if (fTSTimeStep == sampledMinTimeStep)
      {
        /*
         // DEBUG SPECIAL CASE
         if(!reactants)
         {
         G4ExceptionDescription exceptionDescription ;
         exceptionDescription << "No reactants were found by the time stepper.";
         G4Exception("G4Scheduler::ExtractTimeStepperData","ITScheduler008",
         FatalErrorInArgument,exceptionDescription);
         continue ;
         }
         */
        if (bool(reactants))
        {
          // G4cout << "*** (2) G4Scheduler::ExtractTimeStepperData insert
          //          reactants for " << GetIT(track)->GetName() << G4endl;
          // G4cout << "bool(reactants) = " << bool(reactants) << G4endl;
          // G4cout << "trackA : " << GetIT(track)->GetName()
          //        << " ("<< track->GetTrackID() << ")" << G4endl;
          // G4cout << reactants->size() << G4endl;
          // G4cout << GetIT(reactants->operator[](0))->GetName() << G4endl;

         // fReactingTracks.insert(make_pair(track, reactants));
          fReactionSet.AddReactions(fTSTimeStep, track, reactants);
          stepper->ResetReactants();
        }
      }
      else
      {
        if (bool(reactants))
        {
          stepper->ResetReactants();
        }
      }
    }
  }

  MP->CleanProcessor();
}
//_________________________________________________________________________

void G4Scheduler::ComputeInteractionLength()
{
  G4TrackManyList* mainList = fTrackContainer.GetMainList();
  G4TrackManyList::iterator it = mainList ->begin();
  G4TrackManyList::iterator end = mainList ->end();

  fpStepProcessor->SetPreviousStepTime(fPreviousTimeStep);

  for (; it != end; it++)
  {
    G4Track * track = *it;

#ifdef DEBUG
    G4cout << "*CIL* " << GetIT(track)->GetName()
        << " ID: " << track->GetTrackID()
        <<  " at time : " << track->GetGlobalTime()
        << G4endl;
#endif

    fpStepProcessor->DefinePhysicalStepLength(track);

    if (fpStepProcessor->GetTrack()) ExtractILData(fpStepProcessor);
  }
}
//_________________________________________________________________________

void G4Scheduler::ExtractILData(G4ITStepProcessor* SP)
{
  G4Track* track = SP->GetTrack();
  if (!track)
  {
    SP->CleanProcessor();
    return;
  }
  if (track->GetTrackStatus() == fStopAndKill)
  {
    fTrackContainer.GetMainList()->pop(track);
    EndTracking(track);
    //return;
  }

  if (IsInf(SP->GetInteractionTime()))
  {
    // G4cout << "!!!!!!!!!!!! IS INF " << track->GetTrackID() << G4endl;
    SP->CleanProcessor();
    return;
  }
  else if (SP->GetInteractionTime() < fILTimeStep - DBL_EPSILON)
  {
    // G4cout << "!!!!!!!!!!!! TEMPS DIFFERENTS "
    //    << track->GetTrackID() << G4endl;

    ResetLeadingTracks();

    fILTimeStep = SP -> GetInteractionTime();

//    G4cout << "Will set leading step to true for time  :"
//           << SP -> GetInteractionTime() << " against fTimeStep : "
//           << G4BestUnit(fILTimeStep, "Time") << " the trackID is : " << track->GetTrackID()
//           << G4endl;

    GetIT(track)->GetTrackingInfo()->SetLeadingStep(true);
    fLeadingTracks.push_back(track);
  }
  else if(fabs(fILTimeStep - SP -> GetInteractionTime()) < DBL_EPSILON )
  {

    // G4cout << "!!!!!!!!!!!! MEME TEMPS " << track->GetTrackID() << G4endl;
    // G4cout << "Will set leading step to true for time  :"
    //        << SP -> GetInteractionTime() << " against fTimeStep : "
    //        << fTimeStep << " the trackID is : " << track->GetTrackID()<< G4endl;
    GetIT(track)->GetTrackingInfo()->SetLeadingStep(true);
    fLeadingTracks.push_back(track);
  }
  // else
  // {
  //     G4cout << "!!!! Bigger time : " << "currentTime : "<<fILTimeStep
  //  << " proposedTime : " << SP -> GetInteractionTime() << G4endl;
  // }

  SP->CleanProcessor();
}
//_________________________________________________________________________

void G4Scheduler::DoIt()

// Call the process having the min step length or just propagate the track on the given time step

// If the track is "leading the step" (ie one of its process has been selected
// as the one having the minimum time step over all tracks and processes),
// it will undergo its selected processes. Otherwise, it will just propagate the track
// on the given time step.

{

#ifdef G4VERBOSE
  if (fVerbose > 3)
  {
#ifdef USE_COLOR
      G4cout << LIGHT_RED;
#endif
    G4cout << "*** G4Scheduler::DoIt ***" << G4endl;
    G4cout << setw(18) << left << "#Name"
    << setw(15) << "trackID"
    << setw(35) << "Position"
    << setw(25) << "Pre step volume"
    << setw(25) << "Post step volume"
    << setw(22) << "Process"
    << G4endl;
#ifdef USE_COLOR
      G4cout << RESET_COLOR;
#endif
  }
#endif

  G4TrackManyList* mainList = fTrackContainer.GetMainList();
  G4TrackManyList::iterator it = mainList->end();
  it--;
  size_t initialSize = mainList->size();

//  G4cout << "initialSize = " << initialSize << G4endl;

  for(size_t i = 0 ; i < initialSize ; i++)
  {

//    G4cout << "i = " << i << G4endl;

    G4Track* track = *it;
    if (!track)
    {
      G4ExceptionDescription exceptionDescription;
      exceptionDescription << "No track was pop back the main track list.";
      G4Exception("G4Scheduler::DoIt", "ITScheduler009",
          FatalErrorInArgument, exceptionDescription);
    }

#ifdef DEBUG
#ifdef USE_COLOR
      G4cout << LIGHT_RED;
#endif
    G4cout << "*DoIt* " << GetIT(track)->GetName()
            << " ID: " << track->GetTrackID()
            <<  " at time : " << track->GetGlobalTime()
            << G4endl;
#ifdef USE_COLOR
      G4cout << RESET_COLOR;
#endif
#endif
//    G4TrackManyList::iterator next_it (it);
//    next_it--;
//    it = next_it;

    it--;

#ifdef G4VERBOSE
    /////
    // PRE STEP VERBOSE
    if(fVerbose > 3)
    {
      G4String volumeName;

      G4TouchableHandle nextTouchable = track->GetNextTouchableHandle();
      G4VPhysicalVolume* volume (0);

      if(nextTouchable
          && (volume = nextTouchable->GetVolume()))
      {
        volumeName = volume->GetName();

        if(volume->IsParameterised() || volume->IsReplicated())
        {
          volumeName+=" ";
          volumeName+=nextTouchable->GetReplicaNumber();
        }
      }
      else
      {
        volumeName = "OutOfWorld";
      }

      G4cout << setw(18) << left << GetIT(track)->GetName()
      << setw(15) << track->GetTrackID()
      << std::setprecision(3)
      << setw(35) << G4String(G4BestUnit(track->GetPosition(), "Length"))
      << setw(25) << volumeName
      << setw(25) << "---"
      << G4endl;
    }
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    MemStat mem_first, mem_second, mem_diff;
#endif

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    mem_first = MemoryUsage();
#endif

    //	G4cout << " G4Scheduler::DoIt -- fTimeStep="
//    << fTimeStep << G4endl;
    fpStepProcessor->Stepping(track, fTimeStep);

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    MemStat mem_intermediaire = MemoryUsage();
    mem_diff = mem_intermediaire-mem_first;
    G4cout << "\t\t >> || MEM || In DoIT with track "
        << track->GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    ExtractDoItData(fpStepProcessor);

#if defined (DEBUG_MEM) && defined (DEBUG_MEM_DOIT)
    mem_second = MemoryUsage();
    mem_diff = mem_second-mem_first;
    G4cout << "\t >> || MEM || In DoIT with track "
        << track->GetTrackID()
        << ", diff is : " << mem_diff << G4endl;
#endif

#ifdef G4VERBOSE
    /////
    // POST STEP VERBOSE
    if(fVerbose > 3)
    {
      G4cout << setw(18) << left << GetIT(track)->GetName()
      << setw(15) << track->GetTrackID()
      << std::setprecision(3)
      << setw(35) << G4String(G4BestUnit(track->GetPosition(), "Length" ))
      << setw(25) << "---";

      G4TouchableHandle nextTouchable = track->GetNextTouchableHandle();
      G4VPhysicalVolume* volume (0);

      if(nextTouchable
          && (volume = nextTouchable->GetVolume()))
      {
        G4String volumeName = volume->GetName();

        if(volume->IsParameterised() || volume->IsReplicated())
        {
          volumeName+=" ";
          volumeName+=nextTouchable->GetReplicaNumber();
        }

        G4cout << setw(25) << volumeName;
      }
      else
      {
        G4cout << setw(25) << "OutOfWorld";
      }
      if(track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep())
      {
        G4cout << setw(22)
            << track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()
        ->GetProcessName();
      }
      else
      {
        G4cout << "---";
      }
      G4cout << G4endl;

      if(fVerbose > 4)
      {
        const G4TrackVector* secondaries = 0;
        if((secondaries = track->GetStep()->GetSecondary()))
        {
          if(secondaries->empty() == false)
          {
            G4cout << "\t\t ---->";
            for(size_t j = 0; j < secondaries->size(); j++)
            {
              G4cout << GetIT((*secondaries)[j])->GetName()
                  << "("<<(*secondaries)[j]->GetTrackID() << ")"<< " ";
            }
            G4cout << G4endl;
          }
        }
      }
    }
#endif
  }
}
//_________________________________________________________________________

void G4Scheduler::ExtractDoItData(G4ITStepProcessor* SP)
{

  G4Track* track = SP->GetTrack();
  if (!track)
  {
    SP->CleanProcessor();
    return;
  }

  G4TrackStatus status = track->GetTrackStatus();

  switch (status)
  {
    case fAlive:
    case fStopButAlive:
    case fSuspend:
    case fPostponeToNextEvent:
    default:
      PushSecondaries(SP);
      break;

    case fStopAndKill:
      fReactionSet.RemoveReactionSet(track);
      PushSecondaries(SP);
      G4TrackList::Pop(track);
      EndTracking(track);
      break;

    case fKillTrackAndSecondaries:
      fReactionSet.RemoveReactionSet(track);
      G4TrackVector* secondaries = SP->GetSecondaries();
      if (secondaries)
      {
        for (size_t i = 0; i < secondaries->size(); i++)
        {
          delete (*secondaries)[i];
        }
        secondaries->clear();
      }
      G4TrackList::Pop(track);
      EndTracking(track);
      break;
  }
  SP->CleanProcessor();
}
//_________________________________________________________________________

void G4Scheduler::PushSecondaries(G4ITStepProcessor* SP)
{
  G4TrackVector* secondaries = SP->GetSecondaries();
  if (!secondaries || secondaries->empty())
  {
    // DEBUG
    //      G4cout << "NO SECONDARIES !!! " << G4endl;
    return;
  }

  // DEBUG
  //    G4cout << "There are secondaries : "<< secondaries -> size() << G4endl ;

  G4TrackVector::iterator secondaries_i = secondaries->begin();

  for (; secondaries_i != secondaries->end(); secondaries_i++)
  {
    G4Track* secondary = *secondaries_i;
    fTrackContainer._PushTrack(secondary);
  }
}

//_________________________________________________________________________

void G4Scheduler::ComputeTrackReaction()
{
//  if (fReactingTracks.empty())
  if (fReactionSet.Empty())
  {
    return;
  }

  if (fITStepStatus == eCollisionBetweenTracks)
  //        if(fInteractionStep == false)
  {
    fpModelProcessor->FindReaction(&fReactionSet,
                                   fTimeStep,
                                   fPreviousTimeStep, fReachedUserTimeLimit);
    // TODO
    // A ne faire uniquement si le temps choisis est celui calculé par le time stepper
    // Sinon utiliser quelque chose comme : fModelProcessor->FindReaction(&fMainList);

    std::vector<G4ITReactionChange*> * reactionInfo_v = fpModelProcessor
        ->GetReactionInfo();
    std::vector<G4ITReactionChange*>::iterator reactionInfo_i = reactionInfo_v
        ->begin();

    for (; reactionInfo_i != reactionInfo_v->end(); reactionInfo_i++)
    {
      G4ITReactionChange* changes = (*reactionInfo_i);
      G4Track* trackA = const_cast<G4Track*>(changes->GetTrackA());
      G4Track* trackB = const_cast<G4Track*>(changes->GetTrackB());

      if (trackA == 0 || trackB == 0 || trackA->GetTrackStatus() == fStopAndKill
          || trackB->GetTrackStatus() == fStopAndKill) continue;

      G4int nbSecondaries = changes->GetNumberOfSecondaries();

      if (fpUserTimeStepAction)
      {
        const vector<G4Track*>* productsVector = changes->GetfSecondary();
        fpUserTimeStepAction->UserReactionAction(*trackA, *trackB,
                                                 *productsVector);
      }

#ifdef G4VERBOSE
      if (fVerbose)
      {
        G4cout << "At time : " << setw(7) << G4BestUnit(fGlobalTime, "Time")
               << " Reaction : " << GetIT(trackA)->GetName() << " ("
               << trackA->GetTrackID() << ") + " << GetIT(trackB)->GetName()
               << " (" << trackB->GetTrackID() << ") -> ";
      }
#endif

      if (nbSecondaries > 0)
      {
        for (int i = 0; i < nbSecondaries; i++)
        {
#ifdef G4VERBOSE
          if (fVerbose && i != 0) G4cout << " + ";
#endif

          G4Track* secondary = changes->GetSecondary(i);
          fTrackContainer._PushTrack(secondary);
          GetIT(secondary)->SetParentID(trackA->GetTrackID(),
                                        trackB->GetTrackID());

          if (secondary->GetGlobalTime() - fGlobalTime > fTimeTolerance)
          {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription
                << "The time of the secondary should not be bigger than the"
                    " current global time."
                << " This may cause synchronization problem. If the process you"
                    " are using required "
                << "such feature please contact the developpers."
                << G4endl
                << "The global time in the step manager : "
                << G4BestUnit(fGlobalTime,"Time")
                << G4endl
                << "The global time of the track : "
                << G4BestUnit(secondary->GetGlobalTime(),"Time")
                << G4endl;

            G4Exception("G4Scheduler::ComputeInteractionBetweenTracks",
                        "ITScheduler010", FatalErrorInArgument,
                        exceptionDescription);
          }

#ifdef G4VERBOSE
          if (fVerbose) G4cout << GetIT(secondary)->GetName() << " ("
                               << secondary->GetTrackID() << ")";
#endif
        }
      }
      else
      {
#ifdef G4VERBOSE
        if (fVerbose) G4cout << "No product";
#endif
      }
#ifdef G4VERBOSE
      if (fVerbose) G4cout << G4endl;
#endif
      if (trackA->GetTrackID() == 0 || trackB->GetTrackID() == 0)
      {
        G4Track* track = 0;
        if (trackA->GetTrackID() == 0) track = trackA;
        else track = trackB;

        G4ExceptionDescription exceptionDescription;
        exceptionDescription
            << "The problem was found for the reaction between tracks :"
            << trackA->GetParticleDefinition()->GetParticleName() << " ("
            << trackA->GetTrackID() << ") & "
            << trackB->GetParticleDefinition()->GetParticleName() << " ("
            << trackB->GetTrackID() << "). \n";

        if (track->GetStep() == 0)
        {
          exceptionDescription << "Also no step was found"
                               << " ie track->GetStep() == 0 \n";
        }

        exceptionDescription << "Parent ID of trackA : "
                             << trackA->GetParentID() << "\n";
        exceptionDescription << "Parent ID of trackB : "
                             << trackB->GetParentID() << "\n";

        exceptionDescription
            << "The ID of one of the reaction track was not setup.";
        G4Exception("G4Scheduler::ComputeInteractionBetweenTracks",
                    "ITScheduler011", FatalErrorInArgument,
                    exceptionDescription);
      }

      if (changes->WereParentsKilled())
      {
        // DEBUG
        // G4cout << "Erasing tracks : "
        //  << "trackA at time : " << G4BestUnit(trackA->GetGlobalTime() , "Time")
        //  << "\t trackB at time : "<< G4BestUnit(trackB->GetGlobalTime(), "Time")
        //  << "\t GlobalTime : " << G4BestUnit(fGlobalTime, "Time")
        //  << G4endl;

        trackA->SetTrackStatus(fStopAndKill);
        trackB->SetTrackStatus(fStopAndKill);

        G4TrackList::Pop(trackA);
        G4TrackList::Pop(trackB);
        EndTracking(trackA);
        EndTracking(trackB);
      }

      delete changes;
    }

    reactionInfo_v->clear(); // never null because point to a non-pointer member of ModelProcessor
  }
  // DEBUG
  //else
  //{
  //    G4cout << "fInteractionStep == true" << G4endl ;
  //}

//  fReactingTracks.clear();
  fReactionSet.CleanAllReaction();
}

//_________________________________________________________________________

void G4Scheduler::AddTrackID(G4Track* track)
{
  if (fNbTracks == 0) fNbTracks = -1;
  track->SetTrackID(fNbTracks);
  fNbTracks--;
}

//_________________________________________________________________________

void G4Scheduler::PushTrack(G4Track* track)
{
  if (fRunning)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "G4Scheduler::PushTrack : You are trying to push tracks while the "
            "ITScheduler is running";
    G4Exception("G4Scheduler::PushTrack", "ITScheduler012",
                FatalErrorInArgument, exceptionDescription);
  }
  fTrackContainer._PushTrack(track);
}
//_________________________________________________________________________

void G4Scheduler::KillTracks()
{
  fTrackContainer.KillTracks();
}
//_________________________________________________________________________

void G4Scheduler::ResetLeadingTracks()
{
  if (fLeadingTracks.empty() == false)
  {
    std::vector<G4Track*>::iterator fLeadingTracks_i = fLeadingTracks.begin();

    while (fLeadingTracks_i != fLeadingTracks.end())
    {
      G4Track* track = *fLeadingTracks_i;
      if (track)
      {
        G4IT* ITrack = GetIT(*fLeadingTracks_i);
        if (ITrack)
        {
          GetIT(*fLeadingTracks_i)->GetTrackingInfo()->SetLeadingStep(false);
        }
      }

      fLeadingTracks_i++;
      continue;
    }

    fLeadingTracks.clear();
  }
}

//_________________________________________________________________________

void G4Scheduler::EndTracking(G4Track* trackToBeKilled)
{
//  fReactionSet.RemoveReactionSet(trackToBeKilled);
  fpTrackingManager->EndTracking(trackToBeKilled);
  fTrackContainer.PushToKill(trackToBeKilled);
}
//_________________________________________________________________________

void G4Scheduler::EndTracking()
{
  if (fRunning)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "End tracking is called while G4Scheduler is still running."
        << G4endl;

    G4Exception("G4Scheduler::EndTracking", "ITScheduler017",
                FatalErrorInArgument, exceptionDescription);
  }

  fTrackContainer.MergeSecondariesWithMainList();

  if (fTrackContainer.MainListsNOTEmpty())
  {
    G4TrackManyList* mainList = fTrackContainer.GetMainList();
    G4TrackManyList::iterator it = mainList->begin();
    G4TrackManyList::iterator end = mainList->end();
    for (; it != end; it++)
    {
      fpTrackingManager->EndTracking(*it);
    }
  }

  if (fTrackContainer.SecondaryListsNOTEmpty()) // should be empty
  {
    G4TrackManyList* secondaries = fTrackContainer.GetSecondariesList();
    G4TrackManyList::iterator it = secondaries->begin();
    G4TrackManyList::iterator end = secondaries->end();

    for (; it != end; it++)
    {
      fpTrackingManager->EndTracking(*it);
    }
  }
}

//_________________________________________________________________________
void G4Scheduler::SetInteractivity(G4ITTrackingInteractivity* interactivity)
{
  fpTrackingInteractivity = interactivity;
  if (fpTrackingManager) fpTrackingManager->SetInteractivity(
      fpTrackingInteractivity);
}

//_________________________________________________________________________
void G4Scheduler::ForceReinitialization()
{
  fInitialized = false;
  Initialize();
}

//_________________________________________________________________________
G4Scheduler::G4Scheduler(const G4Scheduler&) :
    G4VScheduler(), G4VStateDependent(),
fTrackContainer((G4ITTrackHolder&)*G4ITTrackHolder::Instance())

{
  Create();
}

//_________________________________________________________________________
G4Scheduler& G4Scheduler::operator=(const G4Scheduler& right)
{
  if (this != &right)
  {
    Create();
  }
  return *this;
}

size_t G4Scheduler::GetNTracks()
{
  return fTrackContainer.GetNTracks();
}
//_________________________________________________________________________

void G4Scheduler::GetCollisionType(G4String& interactionType)
{
  switch(fITStepStatus)
  {
    case eInteractionWithMedium:
      interactionType = "eInteractionWithMedium";
      break;
    case eCollisionBetweenTracks:
      interactionType = "eCollisionBetweenTracks";
      break;
    default:
      interactionType = "eCollisionBetweenTracks";
      break;
  }
}

