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

#include <cstdlib>
#include <cstring>
#include <iterator>

//============================================================================//

G4ScoringManager* G4TaskRunManager::masterScM = nullptr;
G4TaskRunManager::masterWorlds_t G4TaskRunManager::masterWorlds =
  G4TaskRunManager::masterWorlds_t();
G4TaskRunManager* G4TaskRunManager::fMasterRM = nullptr;

//============================================================================//

namespace
{
//  G4Mutex cmdHandlingMutex;
  G4Mutex scorerMergerMutex;
  G4Mutex runMergerMutex;
  G4Mutex setUpEventMutex;
}  // namespace

//============================================================================//

G4TaskRunManager* G4TaskRunManager::GetMasterRunManager() { return fMasterRM; }

//============================================================================//

G4RunManagerKernel* G4TaskRunManager::GetMasterRunManagerKernel()
{
  return fMasterRM->kernel;
}

//============================================================================//

G4TaskRunManagerKernel* G4TaskRunManager::GetMTMasterRunManagerKernel()
{
  return fMasterRM->MTkernel;
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
  , workTaskGroupTBB(nullptr)
{
  if(task_queue)
    taskQueue = task_queue;

  // override default of 2 from G4MTRunManager
  nworkers = std::thread::hardware_concurrency();

  if(fMasterRM)
  {
    G4Exception("G4TaskRunManager::G4TaskRunManager", "Run0130", FatalException,
                "Another instance of a G4TaskRunManager already exists.");
  }

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
  G4String _nthread_env = G4GetEnv<G4String>("G4FORCENUMBEROFTHREADS", "max");
  for(auto& itr : _nthread_env)
    itr = tolower(itr);

  if(_nthread_env == "max")
    forcedNwokers = G4Threading::G4GetNumberOfCores();
  else
  {
    std::stringstream ss;
    G4int _nthread_val = -1;
    ss << _nthread_env;
    ss >> _nthread_val;
    if(_nthread_val > 0)
      forcedNwokers = _nthread_val;

    if(forcedNwokers > 0)
    {
      nworkers = forcedNwokers;
      G4cout << "\n### Number of threads is forced to " << forcedNwokers
             << " by Environment variable G4FORCENUMBEROFTHREADS.\n"
             << G4endl;
    }
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
  if(workTaskGroupTBB)
    workTaskGroupTBB->join();
  if(workTaskGroup)
    workTaskGroup->join();
  delete workTaskGroupTBB;
  delete workTaskGroup;
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
    std::stringstream ss;
    ss << "\n### Number of threads is forced to " << forcedNwokers
       << " by G4FORCENUMBEROFTHREADS environment variable. G4TaskRunManager::"
       << __FUNCTION__ << "(" << n << ") ignored ###";

    if(verboseLevel > 1)
    {
      G4ExceptionDescription msg;
      msg << ss.str();
      G4Exception("G4TaskRunManager::SetNumberOfThreads(G4int)", "Run0132",
                  JustWarning, msg);
    }
    else
    {
      G4cout << ss.str() << "\n" << G4endl;
    }
    nworkers = forcedNwokers;
  }
  else
  {
    nworkers = n;
    if(poolInitialized)
    {
      std::stringstream ss;
      ss << "\n### Thread-pool already initialized. Resizing  to " << nworkers
         << "threads ###";
      G4cout << ss.str() << "\n" << G4endl;
      GetThreadPool()->resize(n);
    }
  }
}

//============================================================================//

void G4TaskRunManager::Initialize()
{
  if(!threadPool)
    InitializeThreadPool();

  G4RunManager::Initialize();

  // make sure all worker threads are set up.
  G4RunManager::BeamOn(0);
  G4RunManager::SetRunIDCounter(0);
  // G4UImanager::GetUIpointer()->SetIgnoreCmdNotFound(true);
}

//============================================================================//

void G4TaskRunManager::InitializeThreadPool()
{
  if(poolInitialized && threadPool && (workTaskGroup || workTaskGroupTBB))
  {
    G4Exception("G4TaskRunManager::InitializeThreadPool", "Run1040",
                JustWarning, "Threadpool already initialized. Ignoring...");
    return;
  }

  std::stringstream ss;
  ss.fill('=');
  ss << std::setw(90) << "";
  G4cout << "\n" << ss.str() << G4endl;

  PTL::TaskRunManager::SetVerbose(GetVerbose());
  PTL::TaskRunManager::Initialize(nworkers);

  auto task_join = [](WorkerRunSet& _workers, G4WorkerTaskRunManager** _rm) {
    _workers.insert(_rm);
    return _workers;
  };

  // create the joiners
  if(threadPool->using_tbb())
  {
    G4cout << "G4TaskRunManager :: Using TBB..." << G4endl;
    if(!workTaskGroupTBB)
      workTaskGroupTBB = new RunTaskGroupTBB(task_join, threadPool);
  }
  else
  {
    G4cout << "G4TaskRunManager :: Using G4ThreadPool..." << G4endl;
    if(!workTaskGroup)
      workTaskGroup = new RunTaskGroup(threadPool);
  }

  G4cout << ss.str() << "\n" << G4endl;
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
  G4int grainSize = (eventGrainsize == 0) ? threadPool->size() : eventGrainsize;
  grainSize =
    G4GetEnv<G4int>("G4FORCE_GRAINSIZE", grainSize, "Forcing grainsize...");
  if(grainSize == 0)
    grainSize = 1;

  G4int nEvtsPerTask = (numberOfEventToBeProcessed > grainSize)
                         ? (numberOfEventToBeProcessed / grainSize)
                         : 1;

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

  if(fakeRun)
  {
    std::stringstream msg;
    msg << "--> G4TaskRunManager::ComputeNumberOfTasks() --> " << numberOfTasks
        << " tasks with " << numberOfEventsPerTask << " events/task...";

    std::stringstream ss;
    ss.fill('=');
    ss << std::setw(msg.str().length()) << "";
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
  if(fakeRun)
  {
    {
      std::stringstream msg;
      msg << "--> G4TaskRunManager::CreateAndStartWorkers() --> "
          << "Initializing workers...";

      std::stringstream ss;
      ss.fill('=');
      ss << std::setw(msg.str().length()) << "";
      G4cout << "\n"
             << ss.str() << "\n"
             << msg.str() << "\n"
             << ss.str() << "\n"
             << G4endl;
    }
    if(workTaskGroupTBB)
    {
      G4TaskRunManagerKernel::InitCommandStack() = GetCommandStack();

      // init function which return 1 only once
      auto _init = []() {
        static thread_local G4int _once = 0;
        if(_once++ == 0)
        {
          G4TaskRunManagerKernel::InitializeWorker();
          return 1;
        }
        return 0;
      };

      static std::atomic<size_t> _total_init;
      auto _init_task = [=]() {
        if(G4ThisThread::get_id() == G4MTRunManager::GetMasterThreadId())
          return;
        auto _ret = _init();
        _total_init += _ret;
        if(_ret == 0)
          std::this_thread::sleep_for(std::chrono::milliseconds(100));
      };

      while(_total_init + 1 < threadPool->size())
      {
        std::cout << "repeating init..." << std::endl;
        G4TBBTaskGroup<void> _tbb_init_tg(threadPool);
        auto _n = 2 * threadPool->size();
        while(--_n > 0)
          _tbb_init_tg.exec(_init_task);
        _tbb_init_tg.join();
      }
    }

    if(workTaskGroup)
    {
      threadPool->get_queue()->ExecuteOnAllThreads(
        threadPool, G4TaskRunManagerKernel::InitializeWorker);
    }
  }
  else
  {
    ComputeNumberOfTasks();

    {
      std::stringstream msg;
      msg << "--> G4TaskRunManager::CreateAndStartWorkers() --> "
          << "Creating " << numberOfTasks << " tasks with "
          << numberOfEventsPerTask << " events/task...";

      std::stringstream ss;
      ss.fill('=');
      ss << std::setw(msg.str().length()) << "";
      G4cout << "\n"
             << ss.str() << "\n"
             << msg.str() << "\n"
             << ss.str() << "\n"
             << G4endl;
    }

    G4int remaining = numberOfEventToBeProcessed;

    for(G4int nt = 0; nt < numberOfTasks; ++nt)
    {
      if(remaining == 0)
        break;
      G4int evts = numberOfEventsPerTask;
      if(evts > remaining)
        evts = remaining;
      if(nt + 1 == numberOfTasks)
        evts = remaining;
      AddEventTask(nt, evts);
      remaining -= evts;
    }
  }
}

//============================================================================//

void G4TaskRunManager::AddEventTask(G4int, G4int evts)
{
  TIMEMORY_AUTO_TIMER("");

  // check for TBB first since standard work task group could be valid
  if(workTaskGroupTBB)
    taskManager->exec(*workTaskGroupTBB,
                      G4TaskRunManagerKernel::ExecuteWorkerTask, evts);
  else
    taskManager->exec(*workTaskGroup, G4TaskRunManagerKernel::ExecuteWorkerTask,
                      evts);
}

//============================================================================//

void G4TaskRunManager::RefillSeeds()
{
  TIMEMORY_AUTO_TIMER("");
  G4RNGHelper* helper = G4RNGHelper::GetInstance();
  G4int nFill         = 0;
  switch(SeedOncePerCommunication())
  {
    case 0:
      nFill = numberOfEventToBeProcessed - nSeedsFilled;
      break;
    case 1:
      nFill = numberOfEventsPerTask - nSeedsFilled;
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
  TIMEMORY_AUTO_TIMER("");
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
    if(eventModuloDef > 0)
    {
      eventModulo = eventModuloDef;
      if(eventModulo > numberOfEventsPerTask)
      {
        G4int oldMod = eventModulo;
        eventModulo  = numberOfEventsPerTask;
        if(eventModulo < 1)
          eventModulo = 1;

        G4ExceptionDescription msgd;
        msgd << "Event modulo is reduced to " << eventModulo << " (was "
             << oldMod << ")"
             << " to distribute events to all threads.";
        G4Exception("G4TaskRunManager::InitializeEventLoop()", "Run10035",
                    JustWarning, msgd);
      }
    }
    else
    {
      eventModulo = G4int(std::sqrt(G4double(numberOfEventsPerTask)));
      if(eventModulo < 1)
        eventModulo = 1;
    }

    if(n_event > 0)
    {
      G4bool _overload = InitializeSeeds(n_event);
      G4bool _functor  = false;
      if(!_overload)
        _functor = initSeedsCallback(n_event, nSeedsPerEvent, nSeedsFilled);
      if(_overload == false && _functor == false)
      {
        TIMEMORY_AUTO_TIMER("[reseed]");
        G4RNGHelper* helper = G4RNGHelper::GetInstance();
        switch(SeedOncePerCommunication())
        {
          case 0:
            nSeedsFilled = n_event;
            break;
          case 1:
            nSeedsFilled = numberOfEventsPerTask;
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
            nSeedsFilled             = n_event;
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

  // G4cout << "Preparing command stack..." << G4endl;
  // Prepare UI commands for threads
  PrepareCommandsStack();

  // G4cout << "Create and start workers..." << G4endl;
  // Start worker threads
  CreateAndStartWorkers();

  // G4cout << "Waiting for ready workers..." << G4endl;
  // We need a barrier here. Wait for workers to start event loop.
  // This will return only when all workers have started processing events.
  WaitForReadyWorkers();

  // G4cout << "Event loop initialized..." << G4endl;
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
  size_t nWorlds =
    G4TransportationManager::GetTransportationManager()->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator itrW =
    G4TransportationManager::GetTransportationManager()->GetWorldsIterator();
  for(size_t iWorld = 0; iWorld < nWorlds; ++iWorld)
  {
    addWorld(iWorld, *itrW);
    ++itrW;
  }
}

//============================================================================//

void G4TaskRunManager::MergeScores(const G4ScoringManager* localScoringManager)
{
  TIMEMORY_AUTO_TIMER("");
  G4AutoLock l(&scorerMergerMutex);
  if(masterScM)
    masterScM->Merge(localScoringManager);
}

//============================================================================//

void G4TaskRunManager::MergeRun(const G4Run* localRun)
{
  TIMEMORY_AUTO_TIMER("");
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
    TIMEMORY_AUTO_TIMER("");
    evt->SetEventID(numberOfEventProcessed);
    if(reseedRequired)
    {
      TIMEMORY_AUTO_TIMER("[reseed]");
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
    TIMEMORY_AUTO_TIMER("");
    G4int nev = eventModulo;
    if(numberOfEventProcessed + nev > numberOfEventToBeProcessed)
      nev = numberOfEventToBeProcessed - numberOfEventProcessed;
    evt->SetEventID(numberOfEventProcessed);
    if(reseedRequired)
    {
      TIMEMORY_AUTO_TIMER("[reseed]");
      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      G4int nevRnd        = nev;
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
    numberOfEventProcessed += nev;
    return nev;
  }
  return 0;
}

//============================================================================//

void G4TaskRunManager::TerminateWorkers()
{
  // Force workers to execute (if any) all UI commands left in the stack
  RequestWorkersProcessCommandsStack();
  // Ask workers to exit
  NewActionRequest(WorkerActionRequest::ENDWORKER);

  auto _terminate = [&](WorkerRunSet _workers) {
    for(auto& itr : _workers)
      G4TaskRunManagerKernel::TerminateWorker(itr);
  };

  // Now join threads.
  if(workTaskGroupTBB)
    _terminate(workTaskGroupTBB->join());

  if(workTaskGroup)
  {
    workTaskGroup->join();
    // G4cout << "Terminating workers... (fakeRun: " << std::boolalpha
    //       << fakeRun << ")..." << G4endl;
    if(!fakeRun)
    {
      auto _f = []() { G4TaskRunManagerKernel::TerminateWorker(); };
      threadPool->get_queue()->ExecuteOnAllThreads(threadPool, _f);
    }
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
  auto _terminate = [](WorkerRunSet _workers) {
    for(auto& itr : _workers)
      G4TaskRunManagerKernel::TerminateWorkerRunEventLoop(itr);
  };

  if(workTaskGroupTBB)
    _terminate(workTaskGroupTBB->join());

  if(workTaskGroup)
  {
    workTaskGroup->join();
    if(!fakeRun)
    {
      auto _f = []() { G4TaskRunManagerKernel::TerminateWorkerRunEventLoop(); };
      threadPool->get_queue()->ExecuteOnAllThreads(threadPool, _f);
    }
  }
}

//============================================================================//

void G4TaskRunManager::RequestWorkersProcessCommandsStack()
{
  PrepareCommandsStack();

  auto process_commands_stack = []() {
    G4MTRunManager* mrm        = G4MTRunManager::GetMasterRunManager();
    std::vector<G4String> cmds = mrm->GetCommandStack();
    G4UImanager* uimgr         = G4UImanager::GetUIpointer();  // TLS instance
    for(auto itr = cmds.cbegin(); itr != cmds.cend(); ++itr)
    {
      G4cout << "Applying command \"" << *itr << "\" @ " << __FUNCTION__ << ":"
             << __LINE__ << G4endl;
      uimgr->ApplyCommand(*itr);
    }
    mrm->ThisWorkerProcessCommandsStackDone();
  };

  if(workTaskGroup)
    threadPool->get_queue()->ExecuteOnAllThreads(threadPool,
                                                 process_commands_stack);
}

//============================================================================//

void G4TaskRunManager::ThisWorkerProcessCommandsStackDone() {}

//============================================================================//
