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
// class description:
//   This is a class for run control in GEANT4 for multi-threaded runs
//   It extends G4RunManager re-implementing multi-threaded behavior in
//   key methods. See documentation for G4RunManager
//   Users initializes an instance of this class instead of G4RunManager
//   to start a multi-threaded simulation.

#ifndef G4MTRunManager_h
#define G4MTRunManager_h 1

#include "G4RunManager.hh"
#include "G4Threading.hh"
#include "G4RNGHelper.hh"
#include "G4MTBarrier.hh"
#include <list>
#include <map>

class G4MTRunManagerKernel;
class G4ScoringManager;
class G4UserWorkerInitialization;
class G4UserWorkerThreadInitialization;

//TODO: Split random number storage from this class

class G4MTRunManager : public G4RunManager {
public:
    G4MTRunManager();
    virtual ~G4MTRunManager();
    //New method
    void SetNumberOfThreads( G4int n );
    G4int GetNumberOfThreads() const { return nworkers; }
    void SetPinAffinity(G4int n=1);
    G4int GetPinAffinity() const { return pinAffinity; }
public:

    //Inherited methods to re-implement for MT case
    virtual void Initialize();
    virtual void InitializeEventLoop(G4int n_event, const char* macroFile=0, G4int n_select=-1);

    //The following do not do anything for this runmanager
    virtual void TerminateOneEvent();
    virtual void ProcessOneEvent(G4int i_event);
    ////virtual void TerminateEventLoop();
    virtual void ConstructScoringWorlds();
    virtual void RunTermination();

    //The following method should be invoked by G4WorkerRunManager for each event.
    //False is returned if no more event to be processed.
    // Note: G4Event object must be instantiated by a worker thread. In case no more
    //  event remains to be processed, that worker thread must delete that G4Event
    //  object. If a worker runs with its own random number sequence, the boolean flag
    //  reseedRequired should be set to false. This is *NOT* allowed for the first event.
    virtual G4bool SetUpAnEvent(G4Event*, long& s1, long& s2, long& s3, G4bool reseedRequired=true);
    // Same as above method, but the seeds are set only once over "eventModulo" events.
    // The return value shows the number of events the caller Worker has to process
    // (between 1 and eventModulo depending on number of events yet to be processed).
    // G4Event object has the event ID of the first event of this bunch.
    // If zero is returned no more event needs to be processed, and worker thread 
    // must delete that G4Event.
    virtual G4int SetUpNEvents(G4Event*, G4SeedsQueue* seedsQueue, G4bool reseedRequired=true);
    
    //Method called by Initialize() method
protected:
    //Initialize the seeds list, if derived class does not implement this method
    //A default generation will be used (nevents*2 random seeds)
    //Return true if initialization is done.
    virtual G4bool InitializeSeeds( G4int /*nevts*/) { return false; };
    //Adds one seed to the list of seeds
    virtual void PrepareCommandsStack();
    virtual void StoreRNGStatus(const G4String& filenamePrefix );
    virtual void CreateAndStartWorkers();
    //Creates worker threads and signal to start
public:
    std::vector<G4String> GetCommandStack();
    //This method is invoked just before spawning the threads to
    //collect from UI managere the list of commands that threads
    //will execute.
    size_t GetNumberActiveThreads() const { return threads.size(); }
    //Returns number of currently active threads.
    //This number may be different from the number of threads currently
    //in running state (e.g. the number returned by:
    //G4Threading::GetNumberOfActiveWorkerThreads() method).
private:
    // Number of worker threads. To be set by SetNumberOfThreads() method.
    G4int nworkers;
    // Force to use this number regardless of SetNumberOfThreads() method.
    G4int forcedNwokers;
    // Pin Affinity parameter
    G4int pinAffinity;

    //List of workers (i.e. thread)
    typedef std::list<G4Thread*> G4ThreadsList;
    G4ThreadsList threads;
    //List of workers run managers
    //List of all workers run managers
    std::vector<G4String> uiCmdsForWorkers;
    //List of UI commands for workers.
    CLHEP::HepRandomEngine* masterRNGEngine;
    //Pointer to the mastet thread random engine
protected:
    virtual void WaitForReadyWorkers();
    //Master thread barrier:
    //Call this function to block master thread and
    //wait workers to be ready to process work.
    //This function will return only when all
    //workers are ready to perform event loop.
    virtual void WaitForEndEventLoopWorkers();
    //Master thread barrier:
    //Call this function to block master thread and
    //wait workers have finished current event loop.
    //This function will return only when all
    //workers have finished processing events for this run.
protected:
    G4int numberOfEventToBeProcessed;
    virtual void TerminateWorkers();
    //Empty the workersList

public:
    virtual void ThisWorkerReady();
    //Worker threads barrier:
    //This method should be called by each
    //worker when ready to start thread event-loop
    //This method will return only when all workers
    //are ready.
    //static void ThisWorkerFinishWork();
    //Worker threads barrier:
    //This static method should be called by each
    //worker when finish to process events
    virtual void ThisWorkerEndEventLoop();
    //Worker threads barrier:
    //This method should be called by each
    //worker when worker event loop is terminated.
    typedef std::map<G4int,G4VPhysicalVolume*> masterWorlds_t;
    static G4ScoringManager* GetMasterScoringManager() { return masterScM; }
    static masterWorlds_t& GetMasterWorlds() { return masterWorlds; }
    static void addWorld( G4int counter, G4VPhysicalVolume* w) { masterWorlds.insert( std::make_pair(counter,w) ); }
    const CLHEP::HepRandomEngine* getMasterRandomEngine() const { return masterRNGEngine; }
private:
    //Handling of master thread scoring worlds, access to it is needed by workers
    static G4ScoringManager* masterScM;
    static masterWorlds_t masterWorlds;
    //Singleton implementing master thread behavior
    static G4MTRunManager* fMasterRM;
    G4MTRunManagerKernel* MTkernel;
public: // with description
    static G4MTRunManager* GetMasterRunManager();
    // Returns the singleton instance of the run manager common to all threads implementing 
    // the master behavior
    static G4RunManagerKernel* GetMasterRunManagerKernel();
    static G4MTRunManagerKernel* GetMTMasterRunManagerKernel();
    // Returns the singleton instance of the run manager kernel common to all threads

    virtual void SetUserInitialization(G4VUserPhysicsList* userPL);
    virtual void SetUserInitialization(G4VUserDetectorConstruction* userDC);
    virtual void SetUserInitialization(G4UserWorkerInitialization* userInit);
    virtual void SetUserInitialization(G4UserWorkerThreadInitialization* userInit);
    virtual void SetUserInitialization(G4VUserActionInitialization* userInit);
    virtual void SetUserAction(G4UserRunAction* userAction);
    virtual void SetUserAction(G4VUserPrimaryGeneratorAction* userAction);
    virtual void SetUserAction(G4UserEventAction* userAction);
    virtual void SetUserAction(G4UserStackingAction* userAction);
    virtual void SetUserAction(G4UserTrackingAction* userAction);
    virtual void SetUserAction(G4UserSteppingAction* userAction);

public: 
    // To be invoked solely from G4WorkerRunManager to merge the results
    void MergeScores(const G4ScoringManager* localScoringManager);
    void MergeRun(const G4Run* localRun);
    
public:
    //Handling of more than one run per thread
    enum class WorkerActionRequest {
        UNDEFINED ,
        NEXTITERATION ,    // There is another set of UI commands to be executed
        PROCESSUI,         // Process UI commands w/o a /run/beamOn
        ENDWORKER          // Terminate thread, work finished
    };
    virtual void RequestWorkersProcessCommandsStack();
    //Called to force workers to request and process the UI commands stack
    //This will block untill all workers have processed UI commands
    virtual void ThisWorkerProcessCommandsStackDone();
    //Called by workers to signal to master it has completed processing of
    //UI commands
    virtual WorkerActionRequest ThisWorkerWaitForNextAction();
    //Worker thread barrier
    //This method should be used by workers' run manager to wait,
    //after an event loop for the next action to be performed
    // (for example execute a new run)
    //This returns the action to be performed
protected:
    WorkerActionRequest nextActionRequest;
    virtual void NewActionRequest( WorkerActionRequest newRequest );

protected:
    G4int eventModuloDef;
    G4int eventModulo;
    G4int nSeedsUsed;
    G4int nSeedsFilled;
    G4int nSeedsMax;
    G4int nSeedsPerEvent;
    double* randDbl;

    void RefillSeeds();

public:
    inline void SetEventModulo(G4int i=1) { eventModuloDef = i; }
    inline G4int GetEventModulo() const { return eventModuloDef; }

public:
    virtual void AbortRun(G4bool softAbort=false);
    virtual void AbortEvent();

protected:
    static G4int seedOncePerCommunication;
  // - If it is set to 0 (default), seeds that are centrally managed 
  //   by G4MTRunManager are set for every event of every worker thread. 
  //   This option guarantees event reproducability regardless of number 
  //   of threads. 
  // - If it is set to 1, seeds are set only once for the first 
  //   event of each run of each worker thread. Event reproducability is 
  //   guaranteed only if the same number of worker threads are used. 
  //   On the other hand, this option offers better computing performance 
  //   in particular for applications with relatively small primary 
  //   particle energy and large number of events. 
  // - If it is set to 2, seeds are set only for the first event of 
  //   group of N events. This option is reserved for the future use when 
  //   Geant4 allows number of threads to be dynatically changed during an 
  //   event loop. 
public:
    static G4int SeedOncePerCommunication() { return seedOncePerCommunication; }
    static void SetSeedOncePerCommunication(G4int val) { seedOncePerCommunication = val; }
protected:
    //Barriers: synch points between master and workers
    G4MTBarrier beginOfEventLoopBarrier;
    G4MTBarrier endOfEventLoopBarrier;
    G4MTBarrier nextActionRequestBarrier;
    G4MTBarrier processUIBarrier;
};

#endif //G4MTRunManager_h
