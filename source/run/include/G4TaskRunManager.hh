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
// Author: Jonathan Madsen (May 28st 2020)
//
// class description:
//   This is a class for run control in GEANT4 for multi-threaded runs
//   It extends G4RunManager re-implementing multi-threaded behavior in
//   key methods. See documentation for G4RunManager
//   Users initializes an instance of this class instead of G4RunManager
//   to start a multi-threaded simulation.

#ifndef G4TaskRunManager_hh
#define G4TaskRunManager_hh 1

#include "rundefs.hh"
#include "G4MTBarrier.hh"
#include "G4MTRunManager.hh"
#include "G4RNGHelper.hh"
#include "G4RunManager.hh"
#include "G4TBBTaskGroup.hh"
#include "G4TaskGroup.hh"
#include "G4TaskManager.hh"
#include "G4ThreadPool.hh"
#include "G4Threading.hh"
#include "G4VUserTaskQueue.hh"
#include "G4EnvironmentUtils.hh"
#include "G4Profiler.hh"

#include "PTL/TaskRunManager.hh"

#include <list>
#include <map>

class G4TaskRunManagerKernel;
class G4ScoringManager;
class G4UserTaskInitialization;
class G4UserTaskThreadInitialization;
class G4WorkerTaskRunManager;
class G4RunManagerFactory;

//============================================================================//

class G4TaskRunManager
  : public G4MTRunManager
  , public PTL::TaskRunManager
{
  friend class G4RunManagerFactory;

 public:
  // the profiler aliases are only used when compiled with GEANT4_USE_TIMEMORY
  using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Run>;

 public:
  using InitializeSeedsCallback = std::function<G4bool(G4int, G4int&, G4int&)>;
  using RunTaskGroup            = G4TaskGroup<void>;

 public:
  // Parameters:
  //      taskQueue     : provide a custom task queue
  //      useTBB        : only relevant if GEANT4_USE_TBB defined
  //      evtGrainsize  : the number of events per task
  G4TaskRunManager(G4bool useTBB = G4GetEnv<G4bool>("G4USE_TBB", false));
  G4TaskRunManager(G4VUserTaskQueue* taskQueue,
                   G4bool useTBB      = G4GetEnv<G4bool>("G4USE_TBB", false),
                   G4int evtGrainsize = 0);
  virtual ~G4TaskRunManager();

 public:
  void SetGrainsize(G4int n) { eventGrainsize = n; }
  G4int GetGrainsize() const { return eventGrainsize; }
  inline G4int GetNumberOfTasks() const { return numberOfTasks; }
  inline G4int GetNumberOfEventsPerTask() const
  {
    return numberOfEventsPerTask;
  }

  virtual void SetNumberOfThreads(G4int n) override;
  virtual G4int GetNumberOfThreads() const override
  {
    return PTL::TaskRunManager::GetNumberOfThreads();
  }
  virtual size_t GetNumberActiveThreads() const override
  {
    return PTL::TaskRunManager::GetNumberActiveThreads();
  }
  static G4ThreadId GetMasterThreadId();

 public:
  // Inherited methods to re-implement for MT case
  virtual void Initialize() override;
  virtual void InitializeEventLoop(G4int n_event,
                                   const char* macroFile = nullptr,
                                   G4int n_select        = -1) override;
  virtual void InitializeThreadPool() override;
  G4bool ThreadPoolIsInitialized() const { return poolInitialized; }

  virtual void Initialize(uint64_t nthreads) override
  {
    PTL::TaskRunManager::Initialize(nthreads);
  }

  virtual void TerminateOneEvent() override;
  virtual void ProcessOneEvent(G4int i_event) override;
  virtual void ConstructScoringWorlds() override;
  virtual void RunTermination() override;

  // The following method should be invoked by G4WorkerTaskRunManager for each
  // event. False is returned if no more event to be processed. Note: G4Event
  // object must be instantiated by a worker thread. In case no more
  //  event remains to be processed, that worker thread must delete that G4Event
  //  object. If a worker runs with its own random number sequence, the Boolean
  //  flag reseedRequired should be set to false. This is *NOT* allowed for the
  //  first event.
  virtual G4bool SetUpAnEvent(G4Event*, G4long& s1, G4long& s2, G4long& s3,
                              G4bool reseedRequired = true) override;
  // Same as above method, but the seeds are set only once over "eventModulo"
  // events. The return value shows the number of events the caller Worker has
  // to process (between 1 and eventModulo depending on number of events yet to
  // be processed). G4Event object has the event ID of the first event of this
  // bunch. If zero is returned no more event needs to be processed, and worker
  // thread must delete that G4Event.
  virtual G4int SetUpNEvents(G4Event*, G4SeedsQueue* seedsQueue,
                             G4bool reseedRequired = true) override;

  // Method called by Initialize() method

 protected:
  virtual void ComputeNumberOfTasks();
  // Initialize the seeds list, if derived class does not implement this method
  // A default generation will be used (nevents*2 random seeds)
  // Return true if initialization is done.
  virtual G4bool InitializeSeeds(G4int /*nevts*/) override { return false; }
  virtual void RefillSeeds() override;
  // Adds one seed to the list of seeds
  virtual void StoreRNGStatus(const G4String& filenamePrefix) override;
  virtual void CreateAndStartWorkers() override;
  // Creates worker threads and signal to start

 protected:
  virtual void TerminateWorkers() override;

 public:  // with description
  static G4TaskRunManager* GetMasterRunManager()
  {
    auto* _rm = G4MTRunManager::GetMasterRunManager();
    return dynamic_cast<G4TaskRunManager*>(_rm);
  }
  // Returns the singleton instance of the run manager common to all threads
  // implementing the master behavior
  static G4TaskRunManagerKernel* GetMTMasterRunManagerKernel();
  // Returns the singleton instance of the run manager kernel common to all
  // threads

 public:
  // To be invoked solely from G4WorkerTaskRunManager to merge the results
  void MergeScores(const G4ScoringManager* localScoringManager);
  void MergeRun(const G4Run* localRun);

 public:
  virtual void RequestWorkersProcessCommandsStack() override;
  // Called to force workers to request and process the UI commands stack
  // This will block untill all workers have processed UI commands
  virtual void ThisWorkerProcessCommandsStackDone() override;
  // Called by workers to signal to master it has completed processing of
  // UI commands
  // virtual WorkerActionRequest ThisWorkerWaitForNextAction();
  // Worker thread barrier
  // This method should be used by workers' run manager to wait,
  // after an event loop for the next action to be performed
  // (for example execute a new run)
  // This returns the action to be performed
  virtual void WaitForReadyWorkers() override {}
  virtual void WaitForEndEventLoopWorkers() override;
  virtual void ThisWorkerReady() override {}
  virtual void ThisWorkerEndEventLoop() override {}
  virtual WorkerActionRequest ThisWorkerWaitForNextAction() override
  {
    return WorkerActionRequest::UNDEFINED;
  }

 protected:
  virtual void NewActionRequest(WorkerActionRequest) override {}
  virtual void AddEventTask(G4int);

 public:
  inline void SetInitializeSeedsCallback(InitializeSeedsCallback f)
  {
    initSeedsCallback = f;
  }

 public:
  virtual void AbortRun(G4bool softAbort = false) override;
  virtual void AbortEvent() override;

 private:
  // grainsize
  bool workersStarted                     = false;
  G4int eventGrainsize                    = 0;
  G4int numberOfEventsPerTask             = -1;
  G4int numberOfTasks                     = -1;
  CLHEP::HepRandomEngine* masterRNGEngine = nullptr;
  // Pointer to the master thread random engine
  G4TaskRunManagerKernel* MTkernel = nullptr;

 protected:
  // Barriers: synch points between master and workers
  RunTaskGroup* workTaskGroup = nullptr;

  // aliases to inherited member values
  G4bool& poolInitialized      = PTL::TaskRunManager::m_is_initialized;
  G4ThreadPool*& threadPool    = PTL::TaskRunManager::m_thread_pool;
  G4VUserTaskQueue*& taskQueue = PTL::TaskRunManager::m_task_queue;
  G4TaskManager*& taskManager  = PTL::TaskRunManager::m_task_manager;

  InitializeSeedsCallback initSeedsCallback = [](G4int, G4int&, G4int&) {
    return false;
  };
};

#endif  // G4TaskRunManager_hh
