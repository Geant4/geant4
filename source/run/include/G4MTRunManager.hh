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
// G4MTRunManager
//
// Class description:
//
// This is a class for run control in Geant4 of multi-threaded runs.
// It extends G4RunManager re-implementing multi-threaded behavior in
// key methods (see documentation for G4RunManager).
// Users initialise an instance of this class instead of G4RunManager
// to start a multi-threaded simulation.

// Original authors: X.Dong, A.Dotti - February 2013
// --------------------------------------------------------------------
#ifndef G4MTRunManager_hh
#define G4MTRunManager_hh 1

#include "G4MTBarrier.hh"
#include "G4Profiler.hh"
#include "G4RNGHelper.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"

#include <list>
#include <map>

class G4MTRunManagerKernel;
class G4ScoringManager;
class G4UserWorkerInitialization;
class G4UserWorkerThreadInitialization;
class G4RunManagerFactory;

// TODO: Split random number storage from this class

class G4MTRunManager : public G4RunManager
{
    friend class G4RunManagerFactory;

  public:
    // The profiler aliases are only used when compiled with
    // GEANT4_USE_TIMEMORY.
    using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Run>;

    // Map of defined worlds.
    using masterWorlds_t = std::map<G4int, G4VPhysicalVolume*>;

  public:
    G4MTRunManager();
    ~G4MTRunManager() override;

    void SetNumberOfThreads(G4int n) override;
    G4int GetNumberOfThreads() const override { return nworkers; }
    void SetPinAffinity(G4int n = 1);
    inline G4int GetPinAffinity() const { return pinAffinity; }

    // Inherited methods to re-implement for MT case
    void Initialize() override;
    void InitializeEventLoop(G4int n_event, const char* macroFile = nullptr,
                             G4int n_select = -1) override;
    virtual void InitializeThreadPool() {}

    // The following do not do anything for this runmanager
    void TerminateOneEvent() override;
    void ProcessOneEvent(G4int i_event) override;
    void ConstructScoringWorlds() override;
    void RunTermination() override;

    // The following method should be invoked by G4WorkerRunManager for each
    // event. False is returned if no more event to be processed.
    // Note: G4Event object must be instantiated by a worker thread.
    // In case no more events remain to be processed, that worker thread must
    // delete that G4Event object. If a worker runs with its own random number
    // sequence, the Boolean flag 'reseedRequired' should be set to false.
    // This is *NOT* allowed for the first event.
    virtual G4bool SetUpAnEvent(G4Event*, G4long& s1, G4long& s2, G4long& s3,
                                G4bool reseedRequired = true);

    // Same as above method, but seeds are set only once over "eventModulo"
    // events. The return value shows the number of events the caller Worker
    // has to process (between 1 and eventModulo depending on number of events
    // yet to be processed). G4Event object has the event ID of the first
    // event of this bunch. If zero is returned no more events need to be
    // processed, and worker thread must delete that G4Event.
    // Called by Initialize() method.
    virtual G4int SetUpNEvents(G4Event*, G4SeedsQueue* seedsQueue, G4bool reseedRequired = true);

    // This method is invoked just before spawning the threads to
    // collect from UI manager the list of commands that threads
    // will execute.
    std::vector<G4String> GetCommandStack();

    // Returns number of currently active threads.
    // This number may be different from the number of threads currently
    // in running state, e.g. the number returned by:
    // G4Threading::GetNumberOfActiveWorkerThreads() method.
    virtual size_t GetNumberActiveThreads() const { return threads.size(); }

    static G4ThreadId GetMasterThreadId();

    // Worker threads barrier: this method should be called by each
    // worker when ready to start thread event-loop.
    // This method will return only when all workers are ready.
    virtual void ThisWorkerReady();

    // Worker threads barrier: this method should be called by each
    // worker when worker event loop is terminated.
    virtual void ThisWorkerEndEventLoop();

    static G4ScoringManager* GetMasterScoringManager();
    static masterWorlds_t& GetMasterWorlds();
    static void addWorld(G4int counter, G4VPhysicalVolume* w);

    inline const CLHEP::HepRandomEngine* getMasterRandomEngine() const { return masterRNGEngine; }

    // Returns the singleton instance of the run manager common to all
    // threads implementing the master behavior
    static G4MTRunManager* GetMasterRunManager();

    // Returns the singleton instance of the run manager kernel common to all
    // threads
    static G4RunManagerKernel* GetMasterRunManagerKernel();
    static G4MTRunManagerKernel* GetMTMasterRunManagerKernel();

    void SetUserInitialization(G4VUserPhysicsList* userPL) override;
    void SetUserInitialization(G4VUserDetectorConstruction* userDC) override;
    void SetUserInitialization(G4UserWorkerInitialization* userInit) override;
    void SetUserInitialization(G4UserWorkerThreadInitialization* userInit) override;
    void SetUserInitialization(G4VUserActionInitialization* userInit) override;
    void SetUserAction(G4UserRunAction* userAction) override;
    void SetUserAction(G4VUserPrimaryGeneratorAction* userAction) override;
    void SetUserAction(G4UserEventAction* userAction) override;
    void SetUserAction(G4UserStackingAction* userAction) override;
    void SetUserAction(G4UserTrackingAction* userAction) override;
    void SetUserAction(G4UserSteppingAction* userAction) override;

    // To be invoked solely from G4WorkerRunManager to merge the results
    void MergeScores(const G4ScoringManager* localScoringManager);
    void MergeRun(const G4Run* localRun);

    // Handling of more than one run per thread
    enum class WorkerActionRequest
    {
      UNDEFINED,
      NEXTITERATION,  // There is another set of UI commands to be executed
      PROCESSUI,  // Process UI commands w/o a /run/beamOn
      ENDWORKER  // Terminate thread, work finished
    };

    // Called to force workers to request and process the UI commands stack
    // This will block untill all workers have processed UI commands
    virtual void RequestWorkersProcessCommandsStack();

    // Called by workers to signal to master it has completed processing of
    // UI commands
    virtual void ThisWorkerProcessCommandsStackDone();

    // Worker thread barrier: this method should be used by workers' run
    // manager to wait, after an event loop for the next action to be
    // performed (for example execute a new run).
    // This returns the action to be performed.
    virtual WorkerActionRequest ThisWorkerWaitForNextAction();

    inline void SetEventModulo(G4int i = 1) { eventModuloDef = i; }
    inline G4int GetEventModulo() const { return eventModuloDef; }

    void AbortRun(G4bool softAbort = false) override;
    void AbortEvent() override;

    static G4int SeedOncePerCommunication();
    static void SetSeedOncePerCommunication(G4int val);

  protected:
    // Initialize the seeds list, if derived class does not implement this
    // method, a default generation will be used (nevents*2 random seeds).
    // Return true if initialization is done.
    // Adds one seed to the list of seeds.
    virtual G4bool InitializeSeeds(G4int /*nevts*/) { return false; };

    virtual void PrepareCommandsStack();
    void StoreRNGStatus(const G4String& filenamePrefix) override;
    void rndmSaveThisRun() override;
    void rndmSaveThisEvent() override;

    // Creates worker threads and signal to start
    virtual void CreateAndStartWorkers();

    // Master thread barrier: call this function to block master thread and
    // wait workers to be ready to process work. This function will return
    // only when all workers are ready to perform event loop.
    virtual void WaitForReadyWorkers();

    // Master thread barrier: call this function to block master thread and
    // wait workers have finished current event loop. This function will
    // return only when all workers have finished processing events for
    // this run.
    virtual void WaitForEndEventLoopWorkers();

    // Empty the workersList.
    virtual void TerminateWorkers();

    virtual void NewActionRequest(WorkerActionRequest newRequest);

    virtual void RefillSeeds();

  protected:
    // Number of worker threads. To be set by SetNumberOfThreads() method.
    G4int nworkers = 2;

    // Force to use this number regardless of SetNumberOfThreads() method.
    G4int forcedNwokers = -1;

    G4int numberOfEventToBeProcessed = 0;

    // Handling of master thread scoring worlds, access is needed by workers
    G4MTRUN_DLL static G4ScoringManager* masterScM;
    G4MTRUN_DLL static masterWorlds_t masterWorlds;
    // Singleton implementing master thread behavior
    G4MTRUN_DLL static G4MTRunManager* fMasterRM;

    WorkerActionRequest nextActionRequest = WorkerActionRequest::UNDEFINED;

    G4int eventModuloDef = 0;
    G4int eventModulo = 1;
    G4int nSeedsUsed = 0;
    G4int nSeedsFilled = 0;
    G4int nSeedsMax = 10000;
    G4int nSeedsPerEvent = 2;
    G4double* randDbl = nullptr;

    static G4ThreadId masterThreadId;

    // - If it is set to 0 (default), seeds that are centrally managed
    //   by G4MTRunManager are set for every event of every worker thread.
    //   This option guarantees event reproducibility regardless of number
    //   of threads.
    // - If it is set to 1, seeds are set only once for the first
    //   event of each run of each worker thread. Event reproducibility is
    //   guaranteed only if the same number of worker threads are used.
    //   On the other hand, this option offers better computing performance
    //   in particular for applications with relatively small primary
    //   particle energy and large number of events.
    // - If it is set to 2, seeds are set only for the first event of
    //   group of N events. This option is reserved for the future use when
    //   Geant4 will allow number of threads to be dynamically changed during
    //   an event loop.
    static G4int seedOncePerCommunication;

    // Barriers: synch points between master and workers
    G4MTBarrier beginOfEventLoopBarrier;
    G4MTBarrier endOfEventLoopBarrier;
    G4MTBarrier nextActionRequestBarrier;
    G4MTBarrier processUIBarrier;

  private:
    // List of workers (i.e. thread)
    using G4ThreadsList = std::list<G4Thread*>;

    // Pin Affinity parameter
    G4int pinAffinity = 0;

    // List of workers run managers
    // List of all workers run managers
    G4ThreadsList threads;

    // List of UI commands for workers.
    std::vector<G4String> uiCmdsForWorkers;

    // Pointer to the master thread random engine
    CLHEP::HepRandomEngine* masterRNGEngine = nullptr;

    G4MTRunManagerKernel* MTkernel = nullptr;
};

#endif  // G4MTRunManager_hh
