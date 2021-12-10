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
// G4RunManager
//
// Class description:
//
// This is a class for run control in Geant4
//
// For the sequential mode of Geant4 application, user must provide his
// own classes derived from the following three abstract classes and register
// them to the RunManager:
//   G4VUserDetectorConstruction       - Detector Geometry, Materials
//   G4VUserPhysicsList                - Particle types and Processes
//   G4VUserPrimaryGeneratorAction     - Event Generator selection
//
// In addition to the above mandatory classes, user can easily customise
// the default functionality of a Geant4 simulation by deriving his own
// classes from the following 5 user-action classes:
//   G4UserRunAction                   - Actions for each Run
//   G4UserEventAction                 - Actions for each Event
//   G4UserStackingAction              - Tracks Stacking selection
//   G4UserTrackingAction              - Actions for each Track
//   G4UserSteppingAction              - Actions for each Step
//
// User may use G4VUserActionInitialization class to instantiate any of
// the six user action-classes (1 mandatory + 6 optional).
// In this case, user's concrete G4VUserActionInitialization should be
// defined to the RunManager.
//
// If in multi-threaed mode, a user must provide his own classes derived
// from the following two abstract classes and register them to the
// G4MTRunManager or G4TaskingRunManager:
//   G4VUserDetectorConstruction       - Detector Geometry, Materials
//   G4VUserPhysicsList                - Particle types and Processes
// In addition, user may optionally specify the following:
//   G4UserWorkerInitialization        - Defining thread-local actions
//   G4UserRunAction                   - Actions for entire Run
//
// In multi-threaded mode, the use of G4VUserActionInitialization
// is mandatory.
// In G4VUserActionInitialization, the user has to specify
// G4VUserPrimaryGeneratorAction class. In addition, the default
// functionality of a Geant4 simulation can be customised by making
// user's classes derived from the following 5 user-action classes:
//   G4VUserPrimaryGeneratorAction     - Event Generator selection
//   G4UserRunAction                   - Actions for each tread-local Run
//   G4UserEventAction                 - Actions for each Event
//   G4UserStackingAction              - Tracks Stacking selection
//   G4UserTrackingAction              - Actions for each Track
//   G4UserSteppingAction              - Actions for each Step
//
// G4RunManager MUST be constructed (either explicitly or through
// G4RunManagerFactory) by the user in the main() for enabling sequential
// mode operation of a Geant4 application.
//
// In multi-threaded mode, G4MTRunManager is the dedicated run manager
// which the user MUST construct (either explicitly or through
// G4RunManagerFactory) in the main().
//
// Note: G4WorkerRunManager is the run manager for an individual thread,
// and is instantiated automatically; the user does not need to take care
// of instantiating/deleting it.
// Also, the behavior of the run control can be customised by deriving
// a user class from G4RunManager. In this case, the user should directly
// use the provided protected methods in this class for procedures he/she
// does not want to change.
//
// G4RunManager (or a derived class of it) MUST act as a singleton.
// The user MUST NOT construct more than one such object even if there
// are two different concrete implementations, nor its state can be reset
// to zero, once the object has been created.
//
// G4RunManager controls all of state changes. See G4ApplicationState.hh
// in intercoms category for the meaning of each application state.

// Original author: M.Asai, 1996
// --------------------------------------------------------------------
#ifndef G4RunManager_hh
#define G4RunManager_hh 1

#include <algorithm>
#include <list>

#include "rundefs.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManagerKernel.hh"
#include "G4Profiler.hh"

// userAction classes
class G4VUserDetectorConstruction;
class G4VUserPhysicsList;
class G4UserWorkerInitialization;
class G4UserWorkerThreadInitialization;
class G4VUserActionInitialization;
class G4UserRunAction;
class G4VUserPrimaryGeneratorAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Region;
class G4Timer;
class G4RunMessenger;
class G4DCtable;
class G4Run;
class G4PrimaryTransformer;
class G4RunManagerFactory;

class G4RunManager
{
  friend class G4RunManagerFactory;

  public:

    using ProfilerConfig = G4ProfilerConfig<G4ProfileType::Run>;
      // the profiler aliases are only used when compiled with the
      // GEANT4_USE_TIMEMORY flag enabled.

    static G4RunManager* GetRunManager();
      // Static method which returns the singleton pointer of G4RunManager
      // or its derived class.
      // Note this returns the per-thread singleton in case of a
      // multi-threaded build.

    G4RunManager();
    virtual ~G4RunManager();
      // The constructor and the destructor. The user must construct
      // this class object at the beginning of his/her main() and must
      // delete it at the bottom of the main().

    G4RunManager(const G4RunManager&) = delete;
    G4RunManager& operator=(const G4RunManager&) = delete;
      // Forbidden copy constructor and assignment operator.

    virtual void BeamOn(G4int n_event, const char* macroFile = nullptr,
                        G4int n_select = -1);
      // This method starts an event loop of "n_event" events. The condition of
      // Geant4 is examined before starting the event loop. This method must be
      // invoked at 'Idle' state. The state will be changed to 'GeomClosed'
      // during the event loop and will go back to 'Idle' when the loop is over
      // or aborted.
      // In case a string "macroFile" which represents the name of a macro file
      // is provided, the macro file will be executed AT THE END of each event
      // processing. In case "n_select" is greater than zero, at the end of the
      // first "n_select" events, the macro file is executed.

    virtual void Initialize();
      // This method invokes all the necessary initialisation procedures for an
      // event loop. This method must be invoked at the Geant4 'PreInit' state
      // or 'Idle'. The state will be changed to 'Init' during initialization
      // procedures and then changed to 'Idle'.
      // This method invokes two protected methods, InitializeGeometry() and
      // InitializePhysics().
      // After some event loops, the user can invoke this method once again.
      // It is required if the user changes geometry, physics process, and/or
      // cut-off value. If the user forget the second invocation, the BeamOn()
      // method will invoke this method (Note that this feature is not valid
      // for the first initialization).

    virtual void DefineWorldVolume(G4VPhysicalVolume* worldVol,
                                   G4bool topologyIsChanged = true);
      // This method must be invoked if the geometry setup has been changed
      // between runs. The flag "topologyIsChanged" will specify if the geometry
      // topology is different from the original one used in the previous run;
      // if not, it must be set to false, so that the original optimisation and
      // navigation history are preserved. This method is invoked also at
      // initialisation.

    virtual void AbortRun(G4bool softAbort = false);
      // This method safely aborts the current event loop even if an event is
      // in progress. This method is available for 'GeomClosed' and 'EventProc'
      // Geant4 states. The application state will be changed to 'Idle', so that
      // another event loop can be processed.
      // If the "softAbort" flag is true, the event loop is aborted after
      // processing the current event, while the current event is aborted if the
      // flag is set to false.

    virtual void AbortEvent();
      // This method aborts the currently processing event, remaining events
      // in the current event loop will be processed. This method is available
      // only for 'EventProc' application state.

    virtual void InitializeGeometry();
    virtual void InitializePhysics();
      // These methods are invoked from the Initialize() method for the
      // initializations of geometry and physics processes. The user's concrete
      // G4VUserDetectorConstruction class will be accessed from the method
      // InitializeGeometry() and the G4VUserPhysicsList class will be accessed
      // from the method InitializePhysics().

    virtual G4bool ConfirmBeamOnCondition();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event, const char* macroFile = nullptr,
                             G4int n_select = -1);
    virtual void RunTermination();
      // These four methods are invoked from the BeamOn() method and they're
      // invoked in this order.
      // ConfirmBeamOnCondition() method checks if all the necessary
      // initialisations have been done already. If the condition is not
      // satisfied, false is returned and the following three methods will be
      // skipped.
      // The RunInitialization() method initialises a run. e.g., a G4Run class
      // object is constructed in this method.
      // The DoEventLoop() method controls an event loop. Arguments are the same
      // as for the BeamOn() method.
      // Inside the event loop, the two following methods are invoked at the
      // beginning and at the end of each event.
      // The RunTermination() method terminates a run processing. e.g., a G4Run
      // class object is deleted in this method. If the user adopts ODBMS and
      // wants to store the G4Run object, he/she must override this method.

    virtual void InitializeEventLoop(G4int n_event,
                                     const char* macroFile = nullptr,
                                     G4int n_select = -1);
    virtual void ProcessOneEvent(G4int i_event);
    virtual void TerminateOneEvent();
    virtual void TerminateEventLoop();
      // Granular virtual methods invoked from DoEventLoop().

    virtual G4Event* GenerateEvent(G4int i_event);
    virtual void AnalyzeEvent(G4Event* anEvent);
      // These two methods are invoked from DoEventLoop() at the beginning and
      // at the end of each event processing.
      // GenerateEvent() constructs a G4Event class object and invoke the user's
      // G4VUserPrimaryGeneratorAction concrete class. If the user is adopting
      // an ODBMS system and event objects have been created and stored in the
      // data-base, he/she must override this method.
      // AnalyzeEvent() stores an event to a data-base if a concrete
      // G4VPersistentManager class is defined.

    virtual void ConfigureProfilers(const std::vector<std::string>& args = {});
      // This method configures the global fallback query and label generators.
    void ConfigureProfilers(G4int argc, char** argv);
      // Calls the above virtual method.

    virtual void SetNumberOfThreads(G4int) {}
    virtual G4int GetNumberOfThreads() const { return 1; }
      // Dummy methods to dispatch generic inheritance calls from G4RunManager
      // base class.

    void DumpRegion(const G4String& rname) const;
      // Dump information of a region.
    void DumpRegion(G4Region* region = nullptr) const;
      // Dump information of a region.
      // If the pointer is NULL, all regions are shown.

    void GeometryHasBeenModified(G4bool prop = true);
      // This method must be invoked (or equivalent UI command can be used)
      // in case the user changes his/her detector geometry after Initialize()
      // method has been invoked. Then, at the beginning of the next BeamOn(),
      // all necessary geometry optimisations will be made.
      // The parameter "prop" has to be true if this C++ method is directly
      // invoked.

    void ReinitializeGeometry(G4bool destroyFirst = false, G4bool prop = true);
      // This method must be invoked (or equivalent UI command can be used)
      // in case the user needs his/her detector construction has to be
      // re-invoked. Geometry optimisations will be also done.
      // If the first parameter "destroyFirst" is true, G4SolidStore,
      // G4LogicalVolumeStore and G4PhysicalVolumeStore are cleaned up, and
      // thus all solids, logical volumes and physical volumes previously
      // defined are deleted.
      // The second parameter "prop" has to be true if this C++ method is
      // directly invoked.

    inline void PhysicsHasBeenModified() { kernel->PhysicsHasBeenModified(); }
      // This method must be invoked (or equivalent UI command can be used)
      // in case the user changes his/her physics process(es), e.g. (in)activate
      // some processes. Once this method is invoked, regardless of cuts are
      // changed or not, BuildPhysicsTable() of a PhysicsList is invoked for
      // refreshing all physics tables.

    inline void CutOffHasBeenModified()
    {
      G4cerr << "CutOffHasBeenModified becomes obsolete." << G4endl;
      G4cerr << "It is safe to remove invoking this method." << G4endl;
    }

    void ReOptimizeMotherOf(G4VPhysicalVolume*);
      // This method may be used if the orientation and/or size of a
      // particular physical volume has been modified while the rest of the
      // geometries in the world has not been changed. This avoids the
      // full re-optimisation of the entire geometry tree which is forced
      // if GeometryHasBeenModified() method is invoked.

    void ReOptimize(G4LogicalVolume*);
      // Same as above, but the mother logical volume is specified instead.

    inline void SetGeometryToBeOptimized(G4bool vl)
    {
      if(geometryToBeOptimized != vl)
      {
        geometryToBeOptimized = vl;
        kernel->GeometryHasBeenModified();
        kernel->SetGeometryToBeOptimized(vl);
      }
    }
    inline G4bool GetGeometryToBeOptimized()
    {
      return geometryToBeOptimized;
    }

    void GeometryDirectlyUpdated(G4bool val = true)
    {
      geometryDirectlyUpdated = val;
    }

    static G4bool IfGeometryHasBeenDestroyed();
      // This is used only by workers thread to reset RNG engines from files
      // that are event specific. Not implemented for sequential since run seed
      // defines event seeds.

    virtual void ConstructScoringWorlds();

    virtual void rndmSaveThisRun();
    virtual void rndmSaveThisEvent();
    virtual void RestoreRandomNumberStatus(const G4String& fileN);
    virtual void RestoreRndmEachEvent(G4bool) { /* No effect in SEQ */ }

    virtual void SetUserInitialization(G4VUserDetectorConstruction* userInit);
    virtual void SetUserInitialization(G4VUserPhysicsList* userInit);
    virtual void SetUserInitialization(G4VUserActionInitialization* userInit);
    virtual void SetUserInitialization(G4UserWorkerInitialization* userInit);
    virtual void SetUserInitialization(G4UserWorkerThreadInitialization* userInit);
    virtual void SetUserAction(G4UserRunAction* userAction);
    virtual void SetUserAction(G4VUserPrimaryGeneratorAction* userAction);
    virtual void SetUserAction(G4UserEventAction* userAction);
    virtual void SetUserAction(G4UserStackingAction* userAction);
    virtual void SetUserAction(G4UserTrackingAction* userAction);
    virtual void SetUserAction(G4UserSteppingAction* userAction);
      // Set user-actions and user-initialization to the kernel.
      // Store respective user initialization and action classes.
      // In MT mode, actions are shared among all threads, and should be set
      // in the master thread, while user-actions are thread-private and each      `
      // thread has private instances. Master thread does not have user-actions
      // except for the (optional) run-action.
      // User should instantiate the user-actions in the action-initialization
      // and use that class' setters to set user-actions and *not* directly
      // the methods provided here.
      // Multiple Run, Event, Tracking and Stepping actions are allowed, the
      // multiple instances will be appended to the current configuration.
      // Multiple Stacking and PrimaryGeneration are not allowed.

    inline const G4VUserDetectorConstruction* GetUserDetectorConstruction() const
    {
      return userDetector;
    }
    inline const G4VUserPhysicsList* GetUserPhysicsList() const
    {
      return physicsList;
    }
    inline const G4VUserActionInitialization* GetUserActionInitialization() const
    {
      return userActionInitialization;
    }
    inline G4VUserActionInitialization* GetNonConstUserActionInitialization() const
    {
      return userActionInitialization;
    }
    inline const G4UserWorkerInitialization* GetUserWorkerInitialization() const
    {
      return userWorkerInitialization;
    }
    inline const G4UserWorkerThreadInitialization* GetUserWorkerThreadInitialization() const
    {
      return userWorkerThreadInitialization;
    }
    inline const G4UserRunAction* GetUserRunAction() const
    {
      return userRunAction;
    }
    inline const G4VUserPrimaryGeneratorAction* GetUserPrimaryGeneratorAction() const
    {
      return userPrimaryGeneratorAction;
    }
    inline const G4UserEventAction* GetUserEventAction() const
    {
      return userEventAction;
    }
    inline const G4UserStackingAction* GetUserStackingAction() const
    {
      return userStackingAction;
    }
    inline const G4UserTrackingAction* GetUserTrackingAction() const
    {
      return userTrackingAction;
    }
    inline const G4UserSteppingAction* GetUserSteppingAction() const
    {
      return userSteppingAction;
    }
      // Methods returning respective user initialization and action classes.

    inline void SetNumberOfAdditionalWaitingStacks(G4int iAdd)
      // Set the number of additional (optional) waiting stacks.
      // This method must be invoked at 'PreInit', 'Init' or 'Idle' states.
      // Once the user sets the number of additional waiting stacks,
      // he/she can use the corresponding ENUM in G4ClassificationOfNewTrack.
    {
      eventManager->SetNumberOfAdditionalWaitingStacks(iAdd);
    }

    inline const G4String& GetVersionString() const
    {
      return kernel->GetVersionString();
    }

    inline void SetPrimaryTransformer(G4PrimaryTransformer* pt)
    {
      kernel->SetPrimaryTransformer(pt);
    }

    inline void StoreRandomNumberStatusToG4Event(G4int vl)
      // if vl = 1 : status before primary particle generation is stored
      // if vl = 2 : status before event processing (after primary particle
      // generation) is stored if vl = 3 : both are stored if vl = 0 : none is
      // stored (default).
    {
      storeRandomNumberStatusToG4Event = vl;
      eventManager->StoreRandomNumberStatusToG4Event(vl);
    }

    inline G4int GetFlagRandomNumberStatusToG4Event() const
    {
      return storeRandomNumberStatusToG4Event;
    }

    inline void SetRandomNumberStore(G4bool flag)
    {
      storeRandomNumberStatus = flag;
    }
    inline G4bool GetRandomNumberStore() const
    {
      return storeRandomNumberStatus;
    }
    inline void SetRandomNumberStoreDir(const G4String& dir)
    {
      G4String dirStr = dir;
      if(dirStr.back() != '/')
        dirStr += "/";
    #ifndef WIN32
      G4String shellCmd = "mkdir -p ";
    #else
      std::replace(dirStr.begin(), dirStr.end(), '/', '\\');
      G4String shellCmd = "if not exist " + dirStr + " mkdir ";
    #endif
      shellCmd += dirStr;
      randomNumberStatusDir = dirStr;
      G4int sysret          = system(shellCmd);
      if(sysret != 0)
      {
        G4String errmsg = "\"" + shellCmd
                 + "\" returns non-zero value. Directory creation failed.";
        G4Exception("GrRunManager::SetRandomNumberStoreDir", "Run0071",
                    JustWarning, errmsg);
        G4cerr << " return value = " << sysret << G4endl;
      }
    }
    inline const G4String& GetRandomNumberStoreDir() const
    {
      return randomNumberStatusDir;
    }
    inline const G4String& GetRandomNumberStatusForThisRun() const
    {
      return randomNumberStatusForThisRun;
    }
    inline const G4String& GetRandomNumberStatusForThisEvent() const
    {
      if(storeRandomNumberStatusToG4Event == 0 ||
         storeRandomNumberStatusToG4Event == 2)
      {
        G4Exception("GrRunManager::SetRandomNumberStoreDir", "Run0072",
                    JustWarning,
                    "Random number status is not available for this event.");
      }
      return randomNumberStatusForThisEvent;
    }
    inline void SetRandomNumberStorePerEvent(G4bool flag)
    {
      rngStatusEventsFlag = flag;
    }
    inline G4bool GetRandomNumberStorePerEvent() const
    {
      return rngStatusEventsFlag;
    }

    inline void SetVerboseLevel(G4int vl)
    {
      verboseLevel = vl;
      kernel->SetVerboseLevel(vl);
    }
    inline G4int GetVerboseLevel() const { return verboseLevel; }
    inline G4int GetPrintProgress() { return printModulo; }
    inline void SetPrintProgress(G4int i) { printModulo = i; }

    inline void SetNumberOfEventsToBeStored(G4int val)
      // Sets the number of events to be kept after processing. That is,
      // "val" previous events can be used with the most recent event for
      // digitizing pileup. "val"+1 previous event is deleted.
      // This method must be invoked before starting the event loop.
    {
      n_perviousEventsToBeStored = val;
    }

    inline const G4Run* GetCurrentRun() const { return currentRun; }
    inline G4Run* GetNonConstCurrentRun() const { return currentRun; }
      // Returns the pointer to the current run. This method is available for
      // 'GeomClosed' and 'EventProc' application states.
    inline const G4Event* GetCurrentEvent() const { return currentEvent; }
      // Returns the pointer to the current event. This method is available for
      // 'EventProc' application state.
    inline const G4Event* GetPreviousEvent(G4int i) const
      // Returns the pointer to the "i" previous event. This method is available
      // for 'EventProc' application state. In case the event loop has not yet
      // reached the requested event, null will be returned. To use this method,
      // SetNumberOfEventsToBeStored() method mentioned above must be invoked
      // previously to the event loop.
    {
      if(i >= 1 && i <= n_perviousEventsToBeStored)
      {
        auto itr = previousEvents->cbegin();
        for(G4int j = 1; j < i; ++j)
        {
          ++itr;
        }
        return *itr;
      }
      return nullptr;
    }
    inline void SetRunIDCounter(G4int i) { runIDCounter = i; }
      // Set the run number counter. Initially, the counter is initialized
      // to zero and incremented by one for every BeamOn().

    inline G4int GetNumberOfParallelWorld() const { return nParallelWorlds; }
    inline void SetNumberOfEventsToBeProcessed(G4int val)
    {
      numberOfEventToBeProcessed = val;
    }
    inline G4int GetNumberOfEventsToBeProcessed() const
    {
      return numberOfEventToBeProcessed;
    }
    inline G4int GetNumberOfSelectEvents() const { return n_select_msg; }
    inline const G4String& GetSelectMacro() const { return selectMacro; }
    inline void SetDCtable(G4DCtable* DCtbl) { DCtable = DCtbl; }

    enum RMType
    {
      sequentialRM,
      masterRM,
      workerRM
    };

    inline RMType GetRunManagerType() const { return runManagerType; }

  protected:

    G4RunManager(RMType rmType);
      // This constructor is called in case of multi-threaded build.

    void CleanUpPreviousEvents();
    void CleanUpUnnecessaryEvents(G4int keepNEvents);
    void StackPreviousEvent(G4Event* anEvent);

    virtual void StoreRNGStatus(const G4String& filenamePrefix);

    void UpdateScoring();
    virtual void DeleteUserInitializations();
      // Called by destructor to delete user detector. Note: the user detector
      // is shared among threads, thus this should be re-implemented in derived
      // classes that implement the worker model.

  protected:

    G4RunManagerKernel* kernel = nullptr;
    G4EventManager* eventManager = nullptr;

    G4VUserDetectorConstruction* userDetector = nullptr;
    G4VUserPhysicsList* physicsList = nullptr;
    G4VUserActionInitialization* userActionInitialization = nullptr;
    G4UserWorkerInitialization* userWorkerInitialization = nullptr;
    G4UserWorkerThreadInitialization* userWorkerThreadInitialization = nullptr;
    G4UserRunAction* userRunAction = nullptr;
    G4VUserPrimaryGeneratorAction* userPrimaryGeneratorAction = nullptr;
    G4UserEventAction* userEventAction = nullptr;
    G4UserStackingAction* userStackingAction = nullptr;
    G4UserTrackingAction* userTrackingAction = nullptr;
    G4UserSteppingAction* userSteppingAction = nullptr;

    G4bool geometryInitialized = false;
    G4bool physicsInitialized = false;
    G4bool runAborted = false;
    G4bool initializedAtLeastOnce = false;
    G4bool geometryToBeOptimized = true;

    G4int runIDCounter = 0;
    G4int verboseLevel = 0;
    G4int printModulo = -1;
    G4Timer* timer = nullptr;
    G4DCtable* DCtable = nullptr;

    G4Run* currentRun = nullptr;
    G4Event* currentEvent = nullptr;
    std::list<G4Event*>* previousEvents = nullptr;
    G4int n_perviousEventsToBeStored = 0;
    G4int numberOfEventToBeProcessed = 0;

    G4bool storeRandomNumberStatus = false;
    G4int storeRandomNumberStatusToG4Event = 0;
    G4String randomNumberStatusDir = "./";
    G4String randomNumberStatusForThisRun = "";
    G4String randomNumberStatusForThisEvent = "";
    G4bool rngStatusEventsFlag = false;

    G4VPhysicalVolume* currentWorld = nullptr;

    G4int nParallelWorlds = 0;

    G4String msgText = " ";
    G4int n_select_msg = -1;
    G4int numberOfEventProcessed = 0;
    G4String selectMacro = "";
    G4bool fakeRun = false;
    G4bool isScoreNtupleWriter = false;

    G4bool geometryDirectlyUpdated = false;

    RMType runManagerType;

    G4RUN_DLL static G4bool fGeometryHasBeenDestroyed;
      // This Boolean flag has to be shared by all derived objects.

  private:

    static G4ThreadLocal G4RunManager* fRunManager;
      // Per-thread static instance of the run manager singleton.

    G4RunMessenger* runMessenger = nullptr;

    std::unique_ptr<ProfilerConfig> masterRunProfiler;
};

#endif
