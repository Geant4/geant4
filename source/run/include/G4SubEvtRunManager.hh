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
// G4SubEvtRunManager
//
// Class description:
//
// This is a class for run control in sub-event parallel mode.
// It extends G4TaskRunManager re-implementing multi-threaded behavior in
// key methods (see documentation for G4MTRunManager).
//
// N.B. G4TaskRunManagerKernel is still used by this G4SubEvtRunManager
// without any modification. -- future work for phase-II developments
//
// Original author: M. Asai (JLAB) - 2024
// --------------------------------------------------------------------
#ifndef G4SubEvtRunManager_hh
#define G4SubEvtRunManager_hh 1

#include "G4TaskRunManager.hh"

#include <atomic>
#include <map>

class G4WorkerSubEvtRunManager;

class G4SubEvtRunManager : public G4TaskRunManager
{
    friend class G4RunManagerFactory;

  public:
    static G4SubEvtRunManager* GetMasterRunManager()
    {
      auto* _rm = G4MTRunManager::GetMasterRunManager();
      return dynamic_cast<G4SubEvtRunManager*>(_rm);
    }

    G4SubEvtRunManager(G4bool useTBB = G4GetEnv<G4bool>("G4USE_TBB", false));
    G4SubEvtRunManager(G4VUserTaskQueue* taskQueue,
                     G4bool useTBB = G4GetEnv<G4bool>("G4USE_TBB", false), G4int evtGrainsize = 0);
    ~G4SubEvtRunManager() override;

    // Inherited methods to re-implement for MT case
    void Initialize() override;
    void Initialize(uint64_t nthreads) override { PTL::TaskRunManager::Initialize(nthreads); }
    void RunInitialization() override;
    void InitializeEventLoop(G4int n_event, const char* macroFile = nullptr,
                             G4int n_select = -1) override;
    void TerminateOneEvent() override;
    void ProcessOneEvent(G4int i_event) override;
    void ConstructScoringWorlds() override;
    void RunTermination() override;
    void CleanUpPreviousEvents() override;

    // The following two methods are not used in sub-event parallel mode.
    G4bool SetUpAnEvent(G4Event*, G4long&, G4long&, G4long&, G4bool = true) override
    { return false; }
    G4int SetUpNEvents(G4Event*, G4SeedsQueue*, G4bool = true) override
    { return -1; }

    void RegisterSubEventType(G4int ty, G4int maxEnt) override;

    void RegisterSubEvtWorker(G4WorkerSubEvtRunManager*,G4int);

    // The following method should be invoked by G4WorkerSubEvtRunManager for each
    // sub-event. This method returns a sub-event. This sub-event object must not
    // be deleted by the worker thread.
    // Null pointer is returned if
    //  - run is in progress but no sub-event is available yet (notReady=true), or
    //  - run is over and no more sub-event is available unless new run starts (notReady=false)
    // If a worker runs with its own random number sequence, the Boolean flag 'reseedRequired'
    // should be set to false. This is *NOT* allowed for the first event.
    // G4Event object should be instantiated by G4WorkerSubEvtRunManager based on the
    // input provided in G4SubEvent.
    const G4SubEvent* GetSubEvent(G4int ty, G4bool& notReady,
              G4long& s1, G4long& s2, G4long& s3, G4bool reseedRequired = true) override;

    // Notify the master thread that processing of sub-event "se" is over.
    // Outcome of processed sub-event must be kept in "evt".
    // "se" is deleted by the master thread, while "evt" must be deleted by the worker thread
    // after invoking this method.
    void SubEventFinished(const G4SubEvent* se,const G4Event* evt) override;

    // Merge local scores to the master
    void MergeScores(const G4ScoringManager* localScoringManager) override;

    // Nothing to do for merging run
    void MergeRun(const G4Run*) override
    {}

    // Returns number of currently active threads.
    // This number may be different from the number of threads currently
    // in running state, e.g. the number returned by:
    // G4Threading::GetNumberOfActiveWorkerThreads() method.
    std::size_t GetNumberActiveThreads() const override
    { return threads.size(); }

    // Worker threads barrier: this method should be called by each
    // worker when ready to start thread event-loop.
    // This method will return only when all workers are ready.
    //   void ThisWorkerReady() override;

    // Worker threads barrier: this method should be called by each
    // worker when worker event loop is terminated.
    //   void ThisWorkerEndEventLoop() override;

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

    // Called to force workers to request and process the UI commands stack
    // This will block untill all workers have processed UI commands
    void RequestWorkersProcessCommandsStack() override;

    // Called by workers to signal to master it has completed processing of
    // UI commands
    void ThisWorkerProcessCommandsStackDone() override;

    // Worker thread barrier: this method should be used by workers' run
    // manager to wait, after an event loop for the next action to be
    // performed (for example execute a new run).
    void WaitForReadyWorkers() override {}
    void WaitForEndEventLoopWorkers() override;

    // This returns the action to be performed.
    WorkerActionRequest ThisWorkerWaitForNextAction() override
    {
      return WorkerActionRequest::UNDEFINED;
    }

    void AbortRun(G4bool softAbort = false) override;
    void AbortEvent() override;

  protected:
    void ComputeNumberOfTasks() override
    {
      // This method is not used for sub-event parallel mode
      numberOfTasks = (G4int)threadPool->size();
      numberOfEventsPerTask = -1;
      eventModulo = -1;
    }

    // Initialize the seeds list, if derived class does not implement this
    // method, a default generation will be used (nevents*2 random seeds).
    // Return true if initialization is done.
    // Adds one seed to the list of seeds.
    G4bool InitializeSeeds(G4int /*nevts*/) override { return false; };

    // Adds one seed to the list of seeds
    void RefillSeeds() override;

    //void PrepareCommandsStack() override;

    // Creates worker threads and signal to start
    void CreateAndStartWorkers() override;

    void TerminateWorkers() override;
    void NewActionRequest(WorkerActionRequest) override {}
    void AddEventTask(G4int) override;

    void SetUpSeedsForSubEvent(G4long& s1, G4long& s2, G4long& s3);

    void MergeTrajectories(const G4SubEvent* se,const G4Event* evt) override;
    void UpdateScoringForSubEvent(const G4SubEvent* se,const G4Event* evt) override;

    void CleanUpUnnecessaryEvents(G4int keepNEvents) override;
    void StackPreviousEvent(G4Event* anEvent) override;

  protected:
    std::atomic<G4bool> runInProgress = false;

  private:
    G4bool trajectoriesToBeMerged = false;
    std::map<G4int,G4int> fSubEvtTypeMap;
    std::map<G4WorkerSubEvtRunManager*,G4int> fWorkerMap;

    G4bool CheckSubEvtTypes();

  public:
    void TrajectoriesToBeMerged(G4bool val=true) override
    { trajectoriesToBeMerged = val; }
};

#endif  // G4SubEvtRunManager_hh
