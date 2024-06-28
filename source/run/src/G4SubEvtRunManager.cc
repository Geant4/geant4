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

#include "G4SubEvtRunManager.hh"

#include "G4AutoLock.hh"
#include "G4EnvironmentUtils.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Run.hh"
#include "G4ScoringManager.hh"
#include "G4StateManager.hh"
#include "G4Task.hh"
#include "G4TaskGroup.hh"
#include "G4TaskManager.hh"
#include "G4TaskRunManagerKernel.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4ThreadPool.hh"
#include "G4Threading.hh"
#include "G4Timer.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"
#include "G4UserRunAction.hh"
#include "G4UserTaskInitialization.hh"
#include "G4UserTaskQueue.hh"
#include "G4UserSubEvtThreadInitialization.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UserEventAction.hh"
#include "G4VScoreNtupleWriter.hh"
#include "G4WorkerSubEvtRunManager.hh"
#include "G4WorkerThread.hh"
#include "G4VVisManager.hh"

#include <cstdlib>
#include <cstring>
#include <iterator>
#include <algorithm>

//============================================================================//

namespace
{
G4Mutex scorerMergerMutex;
G4Mutex accessSubEventMutex;
}  // namespace

//============================================================================//

G4SubEvtRunManager::G4SubEvtRunManager(G4VUserTaskQueue* task_queue, G4bool useTBB, G4int grainsize)
 : G4TaskRunManager(task_queue, useTBB, grainsize)
{ runManagerType = subEventMasterRM; }

//============================================================================//

G4SubEvtRunManager::G4SubEvtRunManager(G4bool useTBB)
 : G4SubEvtRunManager(nullptr, useTBB, 0)
{ runManagerType = subEventMasterRM; }

//============================================================================//

G4SubEvtRunManager::~G4SubEvtRunManager()
{
  // relying all the necessary deletion upon the base class
  // G4TaskRunManager::~G4TaskRunManager()
}

//============================================================================//

void G4SubEvtRunManager::Initialize()
{
  G4bool firstTime = (threadPool == nullptr);
  if (firstTime) G4TaskRunManager::InitializeThreadPool();

  G4RunManager::Initialize();

  // make sure all worker threads are set up.
  G4RunManager::BeamOn(0);
  if (firstTime) G4RunManager::SetRunIDCounter(0);
  // G4UImanager::GetUIpointer()->SetIgnoreCmdNotFound(true);
}

//============================================================================//

void G4SubEvtRunManager::RunInitialization()
{
  G4RunManager::RunInitialization();
  if(!fakeRun) runInProgress = true;
}

//============================================================================//

void G4SubEvtRunManager::ProcessOneEvent(G4int i_event)
{
  currentEvent = GenerateEvent(i_event);
  eventManager->ProcessOneEventForSERM(currentEvent);
  
  // Following two lines should not be executed here, as spawned sub-events may 
  // be still being processed. These methods are invoked when all sub-events belongings
  // to this event are processed and scores of these sub-events are marged to the
  // corresponding master event.
  //AnalyzeEvent(currentEvent);
  //UpdateScoring();

  if (i_event < n_select_msg) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
}

//============================================================================//

void G4SubEvtRunManager::TerminateOneEvent()
{
  // We must serialize access here because workers may access Run/Event vectors
  // and we can't override StackPreviousEvent
  G4AutoLock l(&accessSubEventMutex);
  StackPreviousEvent(currentEvent);
  currentEvent = nullptr;
  ++numberOfEventProcessed;
}

//============================================================================//

void G4SubEvtRunManager::StackPreviousEvent(G4Event* anEvent)
{
  if(n_perviousEventsToBeStored>0) {
    G4ExceptionDescription ed;
    ed << "G4RunManager::SetNumberOfEventsToBeStored() is not supported in sub-event parallel mode.\n"
       << "User may still keep events bu G4EventManager::KeepTheCurrentEvent()";
    G4Exception("G4SubEvtRunManager::StackPreviousEvent","SubEvtRM1200",FatalException,ed);
    return;
  }

  if(anEvent->GetNumberOfRemainingSubEvents()>0)
  // sub-events are still under processing. Event is not yet fully completed.
  {
    currentRun->StoreEvent(anEvent);
  }
  else
  // Event is already completed.
  { 
    // making sure this is the first path
    if(!(anEvent->IsEventCompleted()))
    {
      anEvent->EventCompleted();
      if(userEventAction!=nullptr) userEventAction->EndOfEventAction(anEvent);
      auto pVisManager = G4VVisManager::GetConcreteInstance();
      if (pVisManager) pVisManager->EventReadyForVis(anEvent);
      UpdateScoring(anEvent);
      if(!(anEvent->ToBeKept())) 
      { 
        ReportEventDeletion(anEvent);
        delete anEvent;
      }
      else
      { // we keep this event for post-processing (i.e. for vis)
        currentRun->StoreEvent(anEvent);
      }
    } else {
      G4Exception("G4SubEvtRunManager::StackPreviousEvent","SubEvtRM1209",FatalException,"We should not be here!!");
    }
  }

  CleanUpUnnecessaryEvents(0);
}

//============================================================================//

void G4SubEvtRunManager::CleanUpUnnecessaryEvents(G4int keepNEvents)
{
  // Delete events that are no longer necessary for post
  // processing such as visualization.
  // N.B. If ToBeKept() is true, the pointer of this event is
  // kept in G4Run, and deleted along with the deletion of G4Run.

  if(keepNEvents>0) {
    G4ExceptionDescription ed;
    ed << "G4RunManager::SetNumberOfEventsToBeStored() is not supported in sub-event parallel mode.\n"
       << "User may still keep events bu G4EventManager::KeepTheCurrentEvent()";
    G4Exception("G4SubEvtRunManager::CleanUpUnnecessaryEvents","SubEvtRM1201",FatalException,ed);
    return;
  }

  assert(currentRun!=nullptr);

  auto eventVector = currentRun->GetEventVector();
  if(eventVector==nullptr || eventVector->empty()) return;
  auto eItr = eventVector->cbegin();
  while(eItr != eventVector->cend())
  {
    const G4Event* ev = *eItr;
    if(ev!=nullptr)
    {
      if(!(ev->IsEventCompleted()) && ev->GetNumberOfRemainingSubEvents()==0)
      { // This event has been completed since last time we were here
        ev->EventCompleted();
        if(userEventAction!=nullptr) userEventAction->EndOfEventAction(ev);
        auto pVisManager = G4VVisManager::GetConcreteInstance();
        if (pVisManager) pVisManager->EventReadyForVis(ev);
        UpdateScoring(ev);
        if(!(ev->ToBeKept())) 
        {
          ReportEventDeletion(ev);
          delete ev;
          eItr = eventVector->erase(eItr);
        }
        else
        { // we keep this event for post-processing (i.e. for vis)
          eItr++;
        }
      }
      else if(!(ev->ToBeKept()))
      { // post-processing done. we no longer need this event
        ReportEventDeletion(ev);
        delete ev;
        eItr = eventVector->erase(eItr);
      } 
      else
      { // we still need this event
        eItr++;
      }
    }
    else
    { // ev is a null pointer 
      eItr = eventVector->erase(eItr);
    }
  }   
}

//============================================================================//

void G4SubEvtRunManager::CreateAndStartWorkers()
{
  // Now loop on requested number of workers
  // This will also start the workers
  // Currently we do not allow to change the
  // number of threads: threads area created once
  // Instead of pthread based workers, create tbbTask
  static G4bool initializeStarted = false;

  ComputeNumberOfTasks();

  if (fakeRun) {
    if (initializeStarted) {
      auto initCmdStack = GetCommandStack();
      if (!initCmdStack.empty()) {
        threadPool->execute_on_all_threads([initCmdStack]() {
          for (auto& itr : initCmdStack)
            G4UImanager::GetUIpointer()->ApplyCommand(itr);
          G4WorkerTaskRunManager::GetWorkerRunManager()->DoWork();
        });
      }
    }
    else {
      std::stringstream msg;
      msg << "--> G4SubEvtRunManager::CreateAndStartWorkers() --> "
          << "Initializing workers...";

      std::stringstream ss;
      ss.fill('=');
      ss << std::setw((G4int)msg.str().length()) << "";
      G4cout << "\n" << ss.str() << "\n" << msg.str() << "\n" << ss.str() << "\n" << G4endl;

      G4TaskRunManagerKernel::InitCommandStack() = GetCommandStack();
      threadPool->execute_on_all_threads([]() { G4TaskRunManagerKernel::InitializeWorker(); });
    }
    initializeStarted = true;
  }
  else {
    auto initCmdStack = GetCommandStack();
    if (!initCmdStack.empty()) {
      threadPool->execute_on_all_threads([initCmdStack]() {
        for (auto& itr : initCmdStack)
          G4UImanager::GetUIpointer()->ApplyCommand(itr);
      });
    }

    // cleans up a previous run and events in case a thread
    // does not execute any tasks
    threadPool->execute_on_all_threads([]() { G4TaskRunManagerKernel::ExecuteWorkerInit(); });

    {
      std::stringstream msg;
      msg << "--> G4SubEvtRunManager::CreateAndStartWorkers() --> "
          << "Creating " << numberOfTasks << " tasks with " << numberOfEventsPerTask
          << " events/task...";

      std::stringstream ss;
      ss.fill('=');
      ss << std::setw((G4int)msg.str().length()) << "";
      G4cout << "\n" << ss.str() << "\n" << msg.str() << "\n" << ss.str() << "\n" << G4endl;
    }

    /* TODO (PHASE-II): Better calculation of task/event/subevents
       Currently, number of tasks is equal to number of threads
       and each task has a loop that endlessly asks for next sub-event
       until no additional sub-event is available in the master.
       This is not ideal. We should make each task work only for some limited 
       number of sub-events, and create as many number of tasks as needed
       on the fly during the event loop of the master thread., e.g.
    G4int remaining = numberOfEventToBeProcessed;
    for (G4int nt = 0; nt < numberOfTasks + 1; ++nt) {
      if (remaining > 0) AddEventTask(nt);
      remaining -= numberOfEventsPerTask;
    }
    */
    for(G4int nt = 0; nt < numberOfTasks; ++nt)
    { AddEventTask(nt); }
  }
}

//============================================================================//

void G4SubEvtRunManager::AddEventTask(G4int nt)
{
  if (verboseLevel > 1) G4cout << "Adding task " << nt << " to task-group..." << G4endl;
  workTaskGroup->exec([]() { G4TaskRunManagerKernel::ExecuteWorkerTask(); });
}

//============================================================================//

void G4SubEvtRunManager::RefillSeeds()
{
  G4RNGHelper* helper = G4RNGHelper::GetInstance();
  G4int nFill = 0;
  switch (SeedOncePerCommunication()) {
    case 0:
      nFill = numberOfEventToBeProcessed - nSeedsFilled;
      break;
    case 1:
      nFill = numberOfTasks - nSeedsFilled;
      break;
    case 2:
    default:
      nFill = (numberOfEventToBeProcessed - nSeedsFilled * eventModulo) / eventModulo + 1;
  }
  // Generates up to nSeedsMax seed pairs only.
  if (nFill > nSeedsMax) nFill = nSeedsMax;
  masterRNGEngine->flatArray(nSeedsPerEvent * nFill, randDbl);
  helper->Refill(randDbl, nFill);
  nSeedsFilled += nFill;
}

//============================================================================//

void G4SubEvtRunManager::InitializeEventLoop(G4int n_event, const char* macroFile, G4int n_select)
{
  MTkernel->SetUpDecayChannels();
  numberOfEventToBeProcessed = n_event;
  numberOfEventProcessed = 0;

  if (!fakeRun) {
    nSeedsUsed = 0;
    nSeedsFilled = 0;

    if (verboseLevel > 0) timer->Start();

    n_select_msg = n_select;
    if (macroFile != nullptr) {
      if (n_select_msg < 0) n_select_msg = n_event;

      msgText = "/control/execute ";
      msgText += macroFile;
      selectMacro = macroFile;
    }
    else {
      n_select_msg = -1;
      selectMacro = "";
    }

    ComputeNumberOfTasks();

    // initialize seeds
    // If user did not implement InitializeSeeds,
    // use default: nSeedsPerEvent seeds per event

    if (n_event > 0) {
      G4bool _overload = InitializeSeeds(n_event);
      G4bool _functor = false;
      if (!_overload) _functor = initSeedsCallback(n_event, nSeedsPerEvent, nSeedsFilled);
      if (!_overload && !_functor) {
        G4RNGHelper* helper = G4RNGHelper::GetInstance();
        switch (SeedOncePerCommunication()) {
          case 0:
            nSeedsFilled = n_event;
            break;
          case 1:
            nSeedsFilled = numberOfTasks;
            break;
          case 2:
            nSeedsFilled = n_event / eventModulo + 1;
            break;
          default:
            G4ExceptionDescription msgd;
            msgd << "Parameter value <" << SeedOncePerCommunication()
                 << "> of seedOncePerCommunication is invalid. It is reset "
                    "to 0.";
            G4Exception("G4SubEvtRunManager::InitializeEventLoop()", "Run10036", JustWarning, msgd);
            SetSeedOncePerCommunication(0);
            nSeedsFilled = n_event;
        }

        // Generates up to nSeedsMax seed pairs only.
        if (nSeedsFilled > nSeedsMax) nSeedsFilled = nSeedsMax;
        masterRNGEngine->flatArray(nSeedsPerEvent * nSeedsFilled, randDbl);
        helper->Fill(randDbl, nSeedsFilled, n_event, nSeedsPerEvent);
      }
    }
  }

  // Now initialize workers. Check if user defined a WorkerThreadInitialization
  if (userWorkerThreadInitialization == nullptr)
    userWorkerThreadInitialization = new G4UserSubEvtThreadInitialization();

  // Prepare UI commands for threads
  PrepareCommandsStack();

  // Start worker threads
  CreateAndStartWorkers();
}

//============================================================================//

void G4SubEvtRunManager::RunTermination()
{
  // Wait for all worker threads to have finished the run
  // i.e. wait for them to return from RunTermination()
  // This guarantee that userrunaction for workers has been called

  runInProgress = false;

  //TODO (PHASE-II): do we need this???
  workTaskGroup->wait();

  // Wait now for all threads to finish event-loop
  WaitForEndEventLoopWorkers();

  if(currentRun!=nullptr) CleanUpUnnecessaryEvents(0);

  // Now call base-class method
  G4RunManager::TerminateEventLoop();
  G4RunManager::RunTermination();
}

//============================================================================//

void G4SubEvtRunManager::ConstructScoringWorlds()
{
  masterScM = G4ScoringManager::GetScoringManagerIfExist();
  // Call base class stuff...
  G4RunManager::ConstructScoringWorlds();

  masterWorlds.clear();
  auto nWorlds = (G4int)G4TransportationManager::GetTransportationManager()->GetNoWorlds();
  auto itrW = G4TransportationManager::GetTransportationManager()->GetWorldsIterator();
  for (G4int iWorld = 0; iWorld < nWorlds; ++iWorld) {
    addWorld(iWorld, *itrW);
    ++itrW;
  }
}

//============================================================================//

void G4SubEvtRunManager::MergeScores(const G4ScoringManager* localScoringManager)
{
  G4AutoLock l(&scorerMergerMutex);
  if (masterScM != nullptr) masterScM->Merge(localScoringManager);
}

const G4SubEvent* G4SubEvtRunManager::GetSubEvent(G4int ty, G4bool& notReady,
              G4long& s1, G4long& s2, G4long& s3, G4bool reseedRequired)
{
  G4AutoLock l(&accessSubEventMutex);

// This method is invoked from the worker, the ownership of G4SubEvent object
// remains to the master, i.e. will be deleted by the master thread through
// TerminateSubEvent() method.

  if(currentRun==nullptr)
  {
  // Run has not yet started.
    notReady = true;
    return nullptr;
  }

  auto eventVector = currentRun->GetEventVector();
  // RACE HERE: against:
  // 1 G4Run::StoreEvent(G4Event*) G4Run.cc:80
  // 2 G4RunManager::StackPreviousEvent(G4Event*) G4RunManager.cc:572
  for(auto& ev : *eventVector)
  {
  // looping over stored events
    // RACE HERE: against:
    // 1 G4Run::StoreEvent(G4Event*) G4Run.cc:80 
    // 2 G4RunManager::StackPreviousEvent(G4Event*) G4RunManager.cc:572 
    auto se = const_cast<G4Event*>(ev)->PopSubEvent(ty);
    if(se!=nullptr)
    {
    // Sub-event is found in an event that is already finished its event-loop
      notReady = false;
      if(reseedRequired) SetUpSeedsForSubEvent(s1,s2,s3);
      return se;
    }
  }

  auto sep = eventManager->PopSubEvent(ty);
  if(sep!=nullptr)
  {
  // Sub-event is found in an event that is still in the event loop
    notReady = false;
    if(reseedRequired) SetUpSeedsForSubEvent(s1,s2,s3);
    return sep;
  } else {
  // No sub-event available
    // RACE HERE vs line 345
    if(runInProgress)
    {
    // Run is still in progress. Worker should wait until a sub-event is ready
      notReady = true;
    }
    else
    {
    // Run is over. No more sub-event to come unless new run starts.
      notReady = false;
    }
    return nullptr;
  }
}

void G4SubEvtRunManager::SetUpSeedsForSubEvent(G4long& s1, G4long& s2, G4long& s3)
{
  //TODO (PHASE-II): Seeding scheme for sub-event has to be revisited
  G4RNGHelper* helper = G4RNGHelper::GetInstance();
  G4int idx_rndm = nSeedsPerEvent * nSeedsUsed;
  s1 = helper->GetSeed(idx_rndm);
  s2 = helper->GetSeed(idx_rndm + 1);
  if (nSeedsPerEvent == 3) s3 = helper->GetSeed(idx_rndm + 2);
  ++nSeedsUsed;
  if (nSeedsUsed == nSeedsFilled) RefillSeeds();
}


//============================================================================//
//
//G4bool G4SubEvtRunManager::SetUpAnEvent(G4Event* evt, G4long& s1, G4long& s2, G4long& s3,
//                                      G4bool reseedRequired)
//{
//  G4AutoLock l(&setUpEventMutex);
//  if (numberOfEventProcessed < numberOfEventToBeProcessed) {
//    evt->SetEventID(numberOfEventProcessed);
//    if (reseedRequired) {
//      G4RNGHelper* helper = G4RNGHelper::GetInstance();
//      G4int idx_rndm = nSeedsPerEvent * nSeedsUsed;
//      s1 = helper->GetSeed(idx_rndm);
//      s2 = helper->GetSeed(idx_rndm + 1);
//      if (nSeedsPerEvent == 3) s3 = helper->GetSeed(idx_rndm + 2);
//      ++nSeedsUsed;
//      if (nSeedsUsed == nSeedsFilled) RefillSeeds();
//    }
//    numberOfEventProcessed++;
//    return true;
//  }
//  return false;
//}

//============================================================================//
//
//G4int G4SubEvtRunManager::SetUpNEvents(G4Event* evt, G4SeedsQueue* seedsQueue, G4bool reseedRequired)
//{
//  G4AutoLock l(&setUpEventMutex);
//  if (numberOfEventProcessed < numberOfEventToBeProcessed && !runAborted) {
//    G4int nevt = numberOfEventsPerTask;
//    G4int nmod = eventModulo;
//    if (numberOfEventProcessed + nevt > numberOfEventToBeProcessed) {
//      nevt = numberOfEventToBeProcessed - numberOfEventProcessed;
//      nmod = numberOfEventToBeProcessed - numberOfEventProcessed;
//    }
//    evt->SetEventID(numberOfEventProcessed);
//
//    if (reseedRequired) {
//      G4RNGHelper* helper = G4RNGHelper::GetInstance();
//      G4int nevRnd = nmod;
//      if (SeedOncePerCommunication() > 0) nevRnd = 1;
//      for (G4int i = 0; i < nevRnd; ++i) {
//        seedsQueue->push(helper->GetSeed(nSeedsPerEvent * nSeedsUsed));
//        seedsQueue->push(helper->GetSeed(nSeedsPerEvent * nSeedsUsed + 1));
//        if (nSeedsPerEvent == 3) seedsQueue->push(helper->GetSeed(nSeedsPerEvent * nSeedsUsed + 2));
//        nSeedsUsed++;
//        if (nSeedsUsed == nSeedsFilled) RefillSeeds();
//      }
//    }
//    numberOfEventProcessed += nevt;
//    return nevt;
//  }
//  return 0;
//}

//============================================================================//

void G4SubEvtRunManager::SubEventFinished(const G4SubEvent* se,const G4Event* evt)
{
  G4AutoLock l(&accessSubEventMutex);
  UpdateScoringForSubEvent(se,evt);
  eventManager->TerminateSubEvent(se,evt);
}

//============================================================================//

void G4SubEvtRunManager::UpdateScoringForSubEvent(const G4SubEvent* se,const G4Event* evt)
{
  auto masterEvt = se->GetEvent();
  if(masterEvt==nullptr) {
    G4Exception("G4SubEvtRunManager::UpdateScoringForSubEvent()","SERM0001",
      FatalException,"Pointer of master event is null. PANIC!");
  }

  if(userEventAction) { 
    userEventAction->MergeSubEvent(masterEvt,evt);
  }
  //if (isScoreNtupleWriter) {
  //  G4VScoreNtupleWriter::Instance()->Fill(evt->GetHCofThisEvent(),
  //                                         se->GetEvent()->GetEventID());
  //}

  evt->ScoresRecorded();

  G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
  if (ScM == nullptr) return;
  auto nPar = (G4int)ScM->GetNumberOfMesh();
  if (nPar < 1) return;

  if(verboseLevel>2) {
    G4cout << "merging scores of sub-event belonging to event id #" << masterEvt->GetEventID()
           << " --- sub-event has " << evt->GetHCofThisEvent()->GetCapacity()
           << " hits collections" << G4endl;
  }
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if (HCE == nullptr) return;
  auto nColl = (G4int)HCE->GetCapacity();
  for (G4int i = 0; i < nColl; ++i) {
    G4VHitsCollection* HC = HCE->GetHC(i);
    if (HC != nullptr) ScM->Accumulate(HC);
  }
  
}

//============================================================================//

void G4SubEvtRunManager::CleanUpPreviousEvents()
{
  // Delete all events carried over from previous run.
  // This method is invoked at the beginning of the next run
  // or from the destructor of G4RunManager at the very end of
  // the program.
  auto evItr = previousEvents->cbegin();
  while (evItr != previousEvents->cend()) {
    G4Event* evt = *evItr;
    if (evt != nullptr)
    {
      ReportEventDeletion(evt);
      // remove evt from the event vector of G4Run as well
      if(currentRun!=nullptr)
      {
        auto eventVector = currentRun->GetEventVector();
        auto eItr = std::find(eventVector->cbegin(),eventVector->cend(),evt);
        if(eItr != eventVector->cend()) eventVector->erase(eItr);
      }
      delete evt;
    }
    evItr = previousEvents->erase(evItr);
  }
  if(currentRun!=nullptr)
  {
    auto eventVector = currentRun->GetEventVector();
    if(eventVector==nullptr || eventVector->empty()) return;
    auto eItr = eventVector->cbegin();
    while(eItr != eventVector->cend())
    {
      const G4Event* ev = *eItr;
      if(ev!=nullptr)
      {
        ReportEventDeletion(ev);
        delete ev;
      }
      eItr = eventVector->erase(eItr);
    }
  }
}

//============================================================================//

void G4SubEvtRunManager::TerminateWorkers()
{
  // Force workers to execute (if any) all UI commands left in the stack
  RequestWorkersProcessCommandsStack();

  if (workTaskGroup != nullptr) {
    workTaskGroup->join();
    if (!fakeRun)
      threadPool->execute_on_all_threads([]() { G4TaskRunManagerKernel::TerminateWorker(); });
  }
}

//============================================================================//

void G4SubEvtRunManager::AbortRun(G4bool softAbort)
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();
  if (currentState == G4State_GeomClosed || currentState == G4State_EventProc) {
    runAborted = true;
    MTkernel->BroadcastAbortRun(softAbort);
  }
  else {
    G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
  }
}

//============================================================================//

void G4SubEvtRunManager::AbortEvent()
{
  // nothing to do in the master thread
}

//============================================================================//

void G4SubEvtRunManager::WaitForEndEventLoopWorkers()
{
  if (workTaskGroup != nullptr) {
    workTaskGroup->join();
    if (!fakeRun)
      threadPool->execute_on_all_threads(
        []() { G4TaskRunManagerKernel::TerminateWorkerRunEventLoop(); });
  }
}

//============================================================================//

void G4SubEvtRunManager::RequestWorkersProcessCommandsStack()
{
  PrepareCommandsStack();

  auto process_commands_stack = []() {
    G4MTRunManager* mrm = G4MTRunManager::GetMasterRunManager();
    if (mrm != nullptr) {
      auto cmds = mrm->GetCommandStack();
      for (const auto& itr : cmds)
        G4UImanager::GetUIpointer()->ApplyCommand(itr);  // TLS instance
      mrm->ThisWorkerProcessCommandsStackDone();
    }
  };

  if (threadPool != nullptr) threadPool->execute_on_all_threads(process_commands_stack);
}

//============================================================================//

void G4SubEvtRunManager::ThisWorkerProcessCommandsStackDone() {}

//============================================================================//

void G4SubEvtRunManager::SetUserInitialization(G4UserWorkerInitialization* userInit)
{
  userWorkerInitialization = userInit;
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserInitialization(G4UserWorkerThreadInitialization* userInit)
{
  userWorkerThreadInitialization = userInit;
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserInitialization(G4VUserActionInitialization* userInit)
{
  userActionInitialization = userInit;
  userActionInitialization->BuildForMaster();
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserInitialization(G4VUserDetectorConstruction* userDC)
{
  G4RunManager::SetUserInitialization(userDC);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserInitialization(G4VUserPhysicsList* pl)
{
  pl->InitializeWorker();
  G4RunManager::SetUserInitialization(pl);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserAction(G4UserRunAction* userAction)
{
  G4RunManager::SetUserAction(userAction);
  if (userAction != nullptr) userAction->SetMaster(true);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserAction(G4UserEventAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserAction(G4UserStackingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserAction(G4UserTrackingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4SubEvtRunManager::SetUserAction(G4UserSteppingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}


