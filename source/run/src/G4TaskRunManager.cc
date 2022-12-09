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

#include "G4TaskRunManager.hh"
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
#include "G4ThreadPool.hh"
#include "G4Threading.hh"
#include "G4TiMemory.hh"
#include "G4Timer.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"
#include "G4UserRunAction.hh"
#include "G4UserTaskInitialization.hh"
#include "G4UserTaskThreadInitialization.hh"
#include "G4VUserActionInitialization.hh"
#include "G4WorkerTaskRunManager.hh"
#include "G4WorkerThread.hh"
#include "G4UserTaskQueue.hh"
#include "G4TiMemory.hh"
#include "G4ThreadLocalSingleton.hh"

#include <cstdlib>
#include <cstring>
#include <iterator>

//============================================================================//

namespace
{
  G4Mutex scorerMergerMutex;
  G4Mutex runMergerMutex;
  G4Mutex setUpEventMutex;
}  // namespace

//============================================================================//

G4TaskRunManagerKernel* G4TaskRunManager::GetMTMasterRunManagerKernel()
{
  return GetMasterRunManager()->MTkernel;
}

//============================================================================//

G4TaskRunManager::G4TaskRunManager(G4VUserTaskQueue* task_queue, G4bool useTBB,
                                   G4int grainsize)
  : G4MTRunManager()
  , PTL::TaskRunManager(useTBB)
  , eventGrainsize(grainsize)
  , numberOfEventsPerTask(-1)
  , numberOfTasks(-1)
  , masterRNGEngine(nullptr)
  , workTaskGroup(nullptr)
{
  if(task_queue)
    taskQueue = task_queue;

  // override default of 2 from G4MTRunManager
  nworkers  = G4Threading::G4GetNumberOfCores();
  fMasterRM = this;
  MTkernel  = static_cast<G4TaskRunManagerKernel*>(kernel);

  G4int numberOfStaticAllocators = kernel->GetNumberOfStaticAllocators();
  if(numberOfStaticAllocators > 0)
  {
    G4ExceptionDescription msg1;
    msg1 << "There are " << numberOfStaticAllocators
         << " static G4Allocator objects detected.\n"
         << "In multi-threaded mode, all G4Allocator objects must "
         << "be dynamicly instantiated.";
    G4Exception("G4TaskRunManager::G4TaskRunManager", "Run1035", FatalException,
                msg1);
  }

  G4UImanager::GetUIpointer()->SetMasterUIManager(true);
  masterScM = G4ScoringManager::GetScoringManagerIfExist();

  // use default RandomNumberGenerator if created by user, or create default
  masterRNGEngine = G4Random::getTheEngine();

  numberOfEventToBeProcessed = 0;
  randDbl                    = new G4double[nSeedsPerEvent * nSeedsMax];

  //------------------------------------------------------------------------//
  //      handle threading
  //------------------------------------------------------------------------//
  G4String _nthread_env = G4GetEnv<G4String>("G4FORCENUMBEROFTHREADS", "");
  for(auto& itr : _nthread_env)
    itr = (char)std::tolower(itr);

  if(_nthread_env == "max")
    forcedNwokers = G4Threading::G4GetNumberOfCores();
  else if(!_nthread_env.empty())
  {
    std::stringstream ss;
    G4int _nthread_val = -1;
    ss << _nthread_env;
    ss >> _nthread_val;
    if(_nthread_val > 0)
      forcedNwokers = _nthread_val;

    if(forcedNwokers > 0)
      nworkers = forcedNwokers;
  }

  //------------------------------------------------------------------------//
  //      option for forcing TBB
  //------------------------------------------------------------------------//
#ifdef GEANT4_USE_TBB
  G4int _useTBB = G4GetEnv<G4int>("G4FORCE_TBB", (G4int) useTBB);
  if(_useTBB > 0)
    useTBB = true;
#else
  if(useTBB)
  {
    G4ExceptionDescription msg;
    msg << "TBB was requested but Geant4 was not built with TBB support";
    G4Exception("G4TaskRunManager::G4TaskRunManager(...)", "Run0131",
                JustWarning, msg);
  }
  useTBB = false;
#endif

  // handle TBB
  G4ThreadPool::set_use_tbb(useTBB);
}

//============================================================================//

G4TaskRunManager::G4TaskRunManager(G4bool useTBB)
  : G4TaskRunManager(nullptr, useTBB, 0)
{}

//============================================================================//

G4TaskRunManager::~G4TaskRunManager()
{
  // finalize profiler before shutting down the threads
  G4Profiler::Finalize();

  // terminate all the workers
  G4TaskRunManager::TerminateWorkers();

  // trigger all G4AutoDelete instances
  G4ThreadLocalSingleton<void>::Clear();

  // delete the task-group
  delete workTaskGroup;
  workTaskGroup = nullptr;

  // destroy the thread-pool
  if(threadPool)
    threadPool->destroy_threadpool();

  PTL::TaskRunManager::Terminate();
}

//============================================================================//

G4ThreadId G4TaskRunManager::GetMasterThreadId()
{
  return G4MTRunManager::GetMasterThreadId();
}

//============================================================================//

void G4TaskRunManager::StoreRNGStatus(const G4String& fn)
{
  std::ostringstream os;
  os << randomNumberStatusDir << "G4Master_" << fn << ".rndm";
  G4Random::saveEngineStatus(os.str().c_str());
}

//============================================================================//

void G4TaskRunManager::SetNumberOfThreads(G4int n)
{
  if(forcedNwokers > 0)
  {
    if(verboseLevel > 0)
    {
      G4ExceptionDescription msg;
      msg << "\n### Number of threads is forced to " << forcedNwokers
       << " by G4FORCENUMBEROFTHREADS environment variable. G4TaskRunManager::"
       << __FUNCTION__ << "(" << n << ") ignored ###";
      G4Exception("G4TaskRunManager::SetNumberOfThreads(G4int)", "Run0132",
                  JustWarning, msg);
    }
    nworkers = forcedNwokers;
  }
  else
  {
    nworkers = n;
    if(poolInitialized)
    {
      if(verboseLevel > 0)
      {
         std::stringstream ss;
         ss << "\n### Thread-pool already initialized. Resizing  to " << nworkers
            << "threads ###";
         G4cout << ss.str() << "\n" << G4endl;
      }
      GetThreadPool()->resize(n);
    }
  }
}

//============================================================================//

void G4TaskRunManager::Initialize()
{
  G4bool firstTime = (!threadPool);
  if(firstTime)
    InitializeThreadPool();

  G4RunManager::Initialize();

  // make sure all worker threads are set up.
  G4RunManager::BeamOn(0);
  if(firstTime)
    G4RunManager::SetRunIDCounter(0);
  // G4UImanager::GetUIpointer()->SetIgnoreCmdNotFound(true);
}

//============================================================================//

void G4TaskRunManager::InitializeThreadPool()
{
  if(poolInitialized && threadPool && workTaskGroup)
  {
    G4Exception("G4TaskRunManager::InitializeThreadPool", "Run1040",
                JustWarning, "Threadpool already initialized. Ignoring...");
    return;
  }

  PTL::TaskRunManager::SetVerbose(GetVerbose());
  PTL::TaskRunManager::Initialize(nworkers);

  // create the joiners
  if(!workTaskGroup)
  { workTaskGroup = new RunTaskGroup(threadPool); }

  if(verboseLevel > 0)
  {
    std::stringstream ss;
    ss.fill('=');
    ss << std::setw(90) << "";
    G4cout << "\n" << ss.str() << G4endl;

    if(threadPool->is_tbb_threadpool())
    {
      G4cout << "G4TaskRunManager :: Using TBB..." << G4endl;
    }
    else
    {
      G4cout << "G4TaskRunManager :: Using G4ThreadPool..." << G4endl;
    }

    G4cout << ss.str() << "\n" << G4endl;
  }
}

//============================================================================//

void G4TaskRunManager::ProcessOneEvent(G4int)
{
  // Nothing to do
}

//============================================================================//

void G4TaskRunManager::TerminateOneEvent()
{
  // Nothing to do
}

//============================================================================//

void G4TaskRunManager::ComputeNumberOfTasks()
{
  G4int grainSize = (eventGrainsize == 0)
                  ? (G4int)threadPool->size() : eventGrainsize;
  grainSize =
    G4GetEnv<G4int>("G4FORCE_GRAINSIZE", grainSize, "Forcing grainsize...");
  if(grainSize == 0)
    grainSize = 1;

  G4int nEvtsPerTask = (numberOfEventToBeProcessed > grainSize)
                         ? (numberOfEventToBeProcessed / grainSize)
                         : 1;

  if(eventModuloDef > 0)
  {
    eventModulo = eventModuloDef;
  }
  else
  {
    eventModulo = G4int(std::sqrt(G4double(numberOfEventToBeProcessed)));
    if(eventModulo < 1)
      eventModulo = 1;
  }
  if(eventModulo > nEvtsPerTask)
  {
    G4int oldMod = eventModulo;
    eventModulo  = nEvtsPerTask;

    G4ExceptionDescription msgd;
    msgd << "Event modulo is reduced to " << eventModulo << " (was " << oldMod
         << ")"
         << " to distribute events to all threads.";
    G4Exception("G4TaskRunManager::InitializeEventLoop()", "Run10035",
                JustWarning, msgd);
  }
  nEvtsPerTask = eventModulo;

  if(fakeRun)
    nEvtsPerTask = G4GetEnv<G4int>(
      "G4FORCE_EVENTS_PER_TASK", nEvtsPerTask,
      "Forcing number of events per task (overrides grainsize)...");
  else
    nEvtsPerTask = G4GetEnv<G4int>("G4FORCE_EVENTS_PER_TASK", nEvtsPerTask);

  if(nEvtsPerTask < 1)
    nEvtsPerTask = 1;

  numberOfTasks         = numberOfEventToBeProcessed / nEvtsPerTask;
  numberOfEventsPerTask = nEvtsPerTask;
  eventModulo           = numberOfEventsPerTask;

  if(fakeRun && verboseLevel > 1)
  {
    std::stringstream msg;
    msg << "--> G4TaskRunManager::ComputeNumberOfTasks() --> " << numberOfTasks
        << " tasks with " << numberOfEventsPerTask << " events/task...";

    std::stringstream ss;
    ss.fill('=');
    ss << std::setw((G4int)msg.str().length()) << "";
    G4cout << "\n"
           << ss.str() << "\n"
           << msg.str() << "\n"
           << ss.str() << "\n"
           << G4endl;
  }
}

//============================================================================//

void G4TaskRunManager::CreateAndStartWorkers()
{
  // Now loop on requested number of workers
  // This will also start the workers
  // Currently we do not allow to change the
  // number of threads: threads area created once
  // Instead of pthread based workers, create tbbTask
  static bool initializeStarted = false;

  ComputeNumberOfTasks();

  if(fakeRun)
  {
    if(initializeStarted)
    {
      auto initCmdStack = GetCommandStack();
      if(!initCmdStack.empty())
      {
        threadPool->execute_on_all_threads([initCmdStack]() {
          for(auto& itr : initCmdStack)
            G4UImanager::GetUIpointer()->ApplyCommand(itr);
          G4WorkerTaskRunManager::GetWorkerRunManager()->DoWork();
        });
      }
    }
    else
    {
      std::stringstream msg;
      msg << "--> G4TaskRunManager::CreateAndStartWorkers() --> "
          << "Initializing workers...";

      std::stringstream ss;
      ss.fill('=');
      ss << std::setw((G4int)msg.str().length()) << "";
      G4cout << "\n"
             << ss.str() << "\n"
             << msg.str() << "\n"
             << ss.str() << "\n"
             << G4endl;

      G4TaskRunManagerKernel::InitCommandStack() = GetCommandStack();
      threadPool->execute_on_all_threads(
        []() { G4TaskRunManagerKernel::InitializeWorker(); });
    }
    initializeStarted = true;
  }
  else
  {
    auto initCmdStack = GetCommandStack();
    if(!initCmdStack.empty())
    {
      threadPool->execute_on_all_threads([initCmdStack]() {
        for(auto& itr : initCmdStack)
          G4UImanager::GetUIpointer()->ApplyCommand(itr);
      });
    }

    // cleans up a previous run and events in case a thread
    // does not execute any tasks
    threadPool->execute_on_all_threads(
      []() { G4TaskRunManagerKernel::ExecuteWorkerInit(); });

    {
      std::stringstream msg;
      msg << "--> G4TaskRunManager::CreateAndStartWorkers() --> "
          << "Creating " << numberOfTasks << " tasks with "
          << numberOfEventsPerTask << " events/task...";

      std::stringstream ss;
      ss.fill('=');
      ss << std::setw((G4int)msg.str().length()) << "";
      G4cout << "\n"
             << ss.str() << "\n"
             << msg.str() << "\n"
             << ss.str() << "\n"
             << G4endl;
    }

    G4int remaining = numberOfEventToBeProcessed;
    for(G4int nt = 0; nt < numberOfTasks + 1; ++nt)
    {
      if(remaining > 0)
        AddEventTask(nt);
      remaining -= numberOfEventsPerTask;
    }
    workTaskGroup->wait();
  }
}

//============================================================================//

void G4TaskRunManager::AddEventTask(G4int nt)
{
  if(verboseLevel > 1)
    G4cout << "Adding task " << nt << " to task-group..." << G4endl;
  workTaskGroup->exec([]() { G4TaskRunManagerKernel::ExecuteWorkerTask(); });
}

//============================================================================//

void G4TaskRunManager::RefillSeeds()
{
  G4RNGHelper* helper = G4RNGHelper::GetInstance();
  G4int nFill         = 0;
  switch(SeedOncePerCommunication())
  {
    case 0:
      nFill = numberOfEventToBeProcessed - nSeedsFilled;
      break;
    case 1:
      nFill = numberOfTasks - nSeedsFilled;
      break;
    case 2:
    default:
      nFill = (numberOfEventToBeProcessed - nSeedsFilled * eventModulo) /
                eventModulo +
              1;
  }
  // Generates up to nSeedsMax seed pairs only.
  if(nFill > nSeedsMax)
    nFill = nSeedsMax;
  masterRNGEngine->flatArray(nSeedsPerEvent * nFill, randDbl);
  helper->Refill(randDbl, nFill);
  nSeedsFilled += nFill;
}

//============================================================================//

void G4TaskRunManager::InitializeEventLoop(G4int n_event, const char* macroFile,
                                           G4int n_select)
{
  MTkernel->SetUpDecayChannels();
  numberOfEventToBeProcessed = n_event;
  numberOfEventProcessed     = 0;

  if(!fakeRun)
  {
    nSeedsUsed   = 0;
    nSeedsFilled = 0;

    if(verboseLevel > 0)
      timer->Start();

    n_select_msg = n_select;
    if(macroFile != nullptr)
    {
      if(n_select_msg < 0)
        n_select_msg = n_event;

      msgText = "/control/execute ";
      msgText += macroFile;
      selectMacro = macroFile;
    }
    else
    {
      n_select_msg = -1;
      selectMacro  = "";
    }

    ComputeNumberOfTasks();

    // initialize seeds
    // If user did not implement InitializeSeeds,
    // use default: nSeedsPerEvent seeds per event

    if(n_event > 0)
    {
      G4bool _overload = InitializeSeeds(n_event);
      G4bool _functor  = false;
      if(!_overload)
        _functor = initSeedsCallback(n_event, nSeedsPerEvent, nSeedsFilled);
      if(_overload == false && _functor == false)
      {
        G4RNGHelper* helper = G4RNGHelper::GetInstance();
        switch(SeedOncePerCommunication())
        {
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
            G4Exception("G4TaskRunManager::InitializeEventLoop()", "Run10036",
                        JustWarning, msgd);
            SetSeedOncePerCommunication(0);
            nSeedsFilled = n_event;
        }

        // Generates up to nSeedsMax seed pairs only.
        if(nSeedsFilled > nSeedsMax)
          nSeedsFilled = nSeedsMax;
        masterRNGEngine->flatArray(nSeedsPerEvent * nSeedsFilled, randDbl);
        helper->Fill(randDbl, nSeedsFilled, n_event, nSeedsPerEvent);
      }
    }
  }

  // Now initialize workers. Check if user defined a WorkerThreadInitialization
  if(userWorkerThreadInitialization == nullptr)
    userWorkerThreadInitialization = new G4UserTaskThreadInitialization();

  // Prepare UI commands for threads
  PrepareCommandsStack();

  // Start worker threads
  CreateAndStartWorkers();
}

//============================================================================//

void G4TaskRunManager::RunTermination()
{
  // Wait for all worker threads to have finished the run
  // i.e. wait for them to return from RunTermination()
  // This guarantee that userrunaction for workers has been called

  // Wait now for all threads to finish event-loop
  WaitForEndEventLoopWorkers();
  // Now call base-class methof
  G4RunManager::TerminateEventLoop();
  G4RunManager::RunTermination();
}

//============================================================================//

void G4TaskRunManager::ConstructScoringWorlds()
{
  masterScM = G4ScoringManager::GetScoringManagerIfExist();
  // Call base class stuff...
  G4RunManager::ConstructScoringWorlds();

  masterWorlds.clear();
  G4int nWorlds = (G4int)
    G4TransportationManager::GetTransportationManager()->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator itrW =
    G4TransportationManager::GetTransportationManager()->GetWorldsIterator();
  for(G4int iWorld = 0; iWorld < nWorlds; ++iWorld)
  {
    addWorld(iWorld, *itrW);
    ++itrW;
  }
}

//============================================================================//

void G4TaskRunManager::MergeScores(const G4ScoringManager* localScoringManager)
{
  G4AutoLock l(&scorerMergerMutex);
  if(masterScM)
    masterScM->Merge(localScoringManager);
}

//============================================================================//

void G4TaskRunManager::MergeRun(const G4Run* localRun)
{
  G4AutoLock l(&runMergerMutex);
  if(currentRun)
    currentRun->Merge(localRun);
}

//============================================================================//

G4bool G4TaskRunManager::SetUpAnEvent(G4Event* evt, G4long& s1, G4long& s2,
                                      G4long& s3, G4bool reseedRequired)
{
  G4AutoLock l(&setUpEventMutex);
  if(numberOfEventProcessed < numberOfEventToBeProcessed)
  {
    evt->SetEventID(numberOfEventProcessed);
    if(reseedRequired)
    {
      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      G4int idx_rndm      = nSeedsPerEvent * nSeedsUsed;
      s1                  = helper->GetSeed(idx_rndm);
      s2                  = helper->GetSeed(idx_rndm + 1);
      if(nSeedsPerEvent == 3)
        s3 = helper->GetSeed(idx_rndm + 2);
      ++nSeedsUsed;
      if(nSeedsUsed == nSeedsFilled)
        RefillSeeds();
    }
    numberOfEventProcessed++;
    return true;
  }
  return false;
}

//============================================================================//

G4int G4TaskRunManager::SetUpNEvents(G4Event* evt, G4SeedsQueue* seedsQueue,
                                     G4bool reseedRequired)
{
  G4AutoLock l(&setUpEventMutex);
  if(numberOfEventProcessed < numberOfEventToBeProcessed && !runAborted)
  {
    G4int nevt = numberOfEventsPerTask;
    G4int nmod = eventModulo;
    if(numberOfEventProcessed + nevt > numberOfEventToBeProcessed)
    {
      nevt = numberOfEventToBeProcessed - numberOfEventProcessed;
      nmod = numberOfEventToBeProcessed - numberOfEventProcessed;
    }
    evt->SetEventID(numberOfEventProcessed);

    if(reseedRequired)
    {
      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      G4int nevRnd        = nmod;
      if(SeedOncePerCommunication() > 0)
        nevRnd = 1;
      for(G4int i = 0; i < nevRnd; ++i)
      {
        seedsQueue->push(helper->GetSeed(nSeedsPerEvent * nSeedsUsed));
        seedsQueue->push(helper->GetSeed(nSeedsPerEvent * nSeedsUsed + 1));
        if(nSeedsPerEvent == 3)
          seedsQueue->push(helper->GetSeed(nSeedsPerEvent * nSeedsUsed + 2));
        nSeedsUsed++;
        if(nSeedsUsed == nSeedsFilled)
          RefillSeeds();
      }
    }
    numberOfEventProcessed += nevt;
    return nevt;
  }
  return 0;
}

//============================================================================//

void G4TaskRunManager::TerminateWorkers()
{
  // Force workers to execute (if any) all UI commands left in the stack
  RequestWorkersProcessCommandsStack();

  if(workTaskGroup)
  {
    workTaskGroup->join();
    if(!fakeRun)
      threadPool->execute_on_all_threads(
        []() { G4TaskRunManagerKernel::TerminateWorker(); });
  }
}

//============================================================================//

void G4TaskRunManager::AbortRun(G4bool softAbort)
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState =
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState == G4State_GeomClosed || currentState == G4State_EventProc)
  {
    runAborted = true;
    MTkernel->BroadcastAbortRun(softAbort);
  }
  else
  {
    G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
  }
}

//============================================================================//

void G4TaskRunManager::AbortEvent()
{
  // nothing to do in the master thread
}

//============================================================================//

void G4TaskRunManager::WaitForEndEventLoopWorkers()
{
  if(workTaskGroup)
  {
    workTaskGroup->join();
    if(!fakeRun)
      threadPool->execute_on_all_threads(
        []() { G4TaskRunManagerKernel::TerminateWorkerRunEventLoop(); });
  }
}

//============================================================================//

void G4TaskRunManager::RequestWorkersProcessCommandsStack()
{
  PrepareCommandsStack();

  auto process_commands_stack = []() {
    G4MTRunManager* mrm = G4MTRunManager::GetMasterRunManager();
    if(mrm)
    {
      auto cmds = mrm->GetCommandStack();
      for(const auto& itr : cmds)
        G4UImanager::GetUIpointer()->ApplyCommand(itr);  // TLS instance
      mrm->ThisWorkerProcessCommandsStackDone();
    }
  };

  if(threadPool)
    threadPool->execute_on_all_threads(process_commands_stack);
}

//============================================================================//

void G4TaskRunManager::ThisWorkerProcessCommandsStackDone() {}

//============================================================================//
