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
// $Id: G4RunManager.hh 95232 2016-02-01 14:31:22Z gcosmo $
//
// 

// class description:
//
//      This is a class for run control in GEANT4
// 
//     For the sequential mode of Geant4 application,
//     user must provide his own classes derived from the following
//     three abstract classes and register them to the RunManager. 
//        G4VUserDetectorConstruction       - Detector Geometry, Materials
//        G4VUserPhysicsList                - Particle types and Processes 
//        G4VUserPrimaryGeneratorAction     - Event Generator selection
// 
//     In addition to the above mandatory classes, user can easily 
//     customize of the default functionality of GEANT4 simulation
//     by making his own classes derived from the following 5 user
//     action classes. 
//         G4UserRunAction                   - Actions for each Run
//         G4UserEventAction                 - Actions for each Event
//         G4UserStackingAction              - Tracks Stacking selection
//         G4UserTrackingAction              - Actions for each Track
//         G4UserSteppingAction              - Actions for each Step
//
//     User may use G4VUserActionInitialization class to instantiate
//     any of the six user action classes (1 mandatory + 6 optional).
//     In this case, user's concrete G4VUserActionInitialization should
//     be defined to RunManager.
//     
//     For the multi-threaed mode of Geant4 application,
//     user must provide his own classes derived from the following
//     two abstract classes and register them to the MTRunManager.
//        G4VUserDetectorConstruction       - Detector Geometry, Materials
//        G4VUserPhysicsList                - Particle types and Processes
//     In addition, user may optionally specify the following.
//        G4UserWorkerInitialization       - Defining thread-local actions
//        G4UserRunAction                   - Actions for entire Run
//
//     For the multi-threaded mode, use of G4VUserActionInitialization
//     is mandatory.
//     In G4VUserActionInitialization, the user has to specify
//     G4VUserPrimaryGeneratorAction class. In addition user may
//     customize of the default functionality of GEANT4 simulation
//     by making his own classes derived from the following 5 user
//     action classes.
//        G4VUserPrimaryGeneratorAction     - Event Generator selection
//        G4UserRunAction                   - Actions for each tread-local Run
//        G4UserEventAction                 - Actions for each Event
//        G4UserStackingAction              - Tracks Stacking selection
//        G4UserTrackingAction              - Actions for each Track
//        G4UserSteppingAction              - Actions for each Step
//
//     G4RunManager is the only manager class in Geant4 kernel which 
//     the user MUST construct an object by him/herself in the main() 
//     for sequential mode of Geant4 application.
//
//     In the multi-threaded mode, G4MTRunManager is the dedicated
//     run manager which the user MUST construct an object by him/herself
//     in the main().
//
//      Note) G4WorkerRunManager is the run manager for individual
//      thread, and is instantiated automatically, and the user needs
//      not to take care of instantiating/deleting it.
//
//     Also, G4RunManager is the only manager class in Geant4 kernel
//     which the user CAN derive it to costomize the behavior of the
//     run control. For this case, user should use protected methods
//     provided in this class for procedures he/she does not want to
//     change.
//
//     G4RunManager or the derived class of it MUST be a singleton.
//     The user MUST NOT construct more than one object even if there
//     are two different concrete implementations.
//
//     G4RunManager controls all of state changes. See G4ApplicationState.hh
//     in intercoms category for the meanings of each state.
//

#ifndef G4RunManager_h
#define G4RunManager_h 1

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

#include "G4RunManagerKernel.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "globals.hh"
#include <list>
#include <algorithm>

class G4RunManager
{
  public: // with description
    static G4RunManager* GetRunManager();
    //  Static method which returns the singleton pointer of G4RunManager or
    // its derived class.
    // Note this returns the per-thread singleton in case of multi-threaded
    // build

  private:
    static G4ThreadLocal G4RunManager* fRunManager;
    //Per-thread static instance of the run manager singleton
 
public: // with description
    G4RunManager();

    virtual ~G4RunManager();
    //  The constructor and the destructor. The user must construct this class
    // object at the beginning of his/her main() and must delete it at the 
    // bottom of the main().

  public: // with description
    virtual void BeamOn(G4int n_event,const char* macroFile=0,G4int n_select=-1);
    //  This method starts an event loop of "n_event" events. The condition of Geant4
    // is examined before starting the event loop. This method must be invoked at
    // Idle state. The state will be changed to GeomClosed during the event loop and
    // will go back to Idle when the loop is over or aborted.
    //  In case a string "macroFile" which represents the name of a macro file is given,
    // this macro file will be executed AT THE END of each event processing. In case
    // "n_select" is greater than zero, at the end of first "n_select" events the macro
    // file is executed.
    virtual void Initialize();
    //  This method invokes all the necessary initialization procedures for an event
    // loop. This method must be invoked at the Geant4 state of PreInit or Idle. The
    // state will be changed to Init during the initialization procedures and then
    // changed to Idle.
    //  This method invokes two protected methods, InitializeGeometry() and
    // InitializePhysics().
    //  After some event loops, the user can invoke this method once again. It is
    // required if the user changes geometry, physics process, and/or cut off value.
    // If the user forget the second invokation, G4RunManager will invoke BeamOn()
    // method will invoke this method. (Note that this feature is not valid for the
    // first initialization.)
    virtual void DefineWorldVolume(G4VPhysicalVolume * worldVol,
                                   G4bool topologyIsChanged=true);
    //  This method must be invoked if the geometry setup has been changed between
    // runs. The flag 'topologyIsChanged' will specify if the geometry topology is
    // different from the original one used in the previous run; if not, it must be
    // set to false, so that the original optimisation and navigation history is
    // preserved. This method is invoked also at initialisation.
    //////////////////////////////////////////////////////virtual void ResetNavigator() const;
    //  Resets the state of the navigator for tracking; needed for geometry updates.
    // It forces the optimisation and navigation history to be reset.
    virtual void AbortRun(G4bool softAbort=false);
    //  This method safely aborts the current event loop even if an event is in progress.
    // This method is available for Geant4 states of GeomClosed and EventProc. The state
    // will be changed to Idle, so that another event loop can be done.
    //  If softAbort is true, the event loop is aborted after processing the current
    // event, while the current event is aborted if it is false.
    virtual void AbortEvent();
    //  This method aborts the currently processing event, remaining events in the
    // current event loop will be processed. This method is available only for
    // EventProc state.
  public: // with description

    virtual void InitializeGeometry();
    virtual void InitializePhysics();
    //  These protected methods are invoked from Initialize() method for the 
    // initializations of geometry and physics processes. The user's concrete
    // G4VUserDetectorConstruction class will be accessed from InitializeGeometry() and
    // G4VUserPhysicsList class will be accessed from InitializePhysics().

    virtual G4bool ConfirmBeamOnCondition();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event,const char* macroFile=0,G4int n_select=-1);
    virtual void RunTermination();
    //  These four protected methods are invoked from BeamOn() method. These four methods
    // are invoked in this order.
    //  ConfirmBeamOnCondition() method checks if all the necessary initializations have
    // already done. If the condition is not satisfied, false is returned and the follwing
    // three methods will be skipped.
    //  RunInitialization() method initializes a run. For example, a G4Run class object 
    // is constructed in this method.
    //  DoEventLoop() method control an event loop. Arguments are same as BeamOn() method.
    // Inide the event loop, two following protected methods are invoked at the begining
    // and the end of each event. 
    //  RunTermination() method terminates a run processing. For example, a G4Run class
    // object is deleted in this class. If the user uses ODBMS and wants to store the
    // G4Run class object, he/she must override this method.

    virtual void InitializeEventLoop(G4int n_event,const char* macroFile=0,G4int n_select=-1);
    virtual void ProcessOneEvent(G4int i_event);
    virtual void TerminateOneEvent();
    virtual void TerminateEventLoop();
    //  Granular virtual methods invoked from DoEventLoop() method.

    ///////////////////////////////////////////////////////////virtual void BuildPhysicsTables();
    //  This method is invoked from RunInitialization() to create physics tables.

    virtual G4Event* GenerateEvent(G4int i_event);
    virtual void AnalyzeEvent(G4Event* anEvent);
    //  These two protected methods are invoked from DoEventLoop() method at the begining
    // and the end of each event processing.
    //  GenerateEvent() method constructs a G4Event class object and invoke the user's
    // G4VUserPrimaryGeneratorAction concrete class. If the user is using ODBMS and event 
    // objects have been created and stored in the data base, he/she must override this
    // method.
    //  AnalyzeEvent() stores an event to a data base if a concrete G4VPersistentManager
    // class is defined.

  public: // with description
    //////////////////////////////////////////////////////void UpdateRegion();
    // Update region list. 
    // This method is mandatory before invoking following two dump methods.
    // At RunInitialization(), this method is automatically invoked, and thus
    // the user needs not invoke.

    void DumpRegion(const G4String& rname) const;
    // Dump information of a region.

    void DumpRegion(G4Region* region=0) const;
    // Dump information of a region.
    // If the pointer is NULL, all regions are shown.

  protected:
    void CleanUpPreviousEvents();
    void CleanUpUnnecessaryEvents(G4int keepNEvents);
    void StackPreviousEvent(G4Event* anEvent);

  public:
    enum RMType { sequentialRM, masterRM, workerRM };
  protected:
    //This constructor is called in case of Geant4 Multi-threaded build
    G4RunManager( RMType rmType );

  protected:
    G4RunManagerKernel * kernel;
    G4EventManager * eventManager;

    G4VUserDetectorConstruction * userDetector;
    G4VUserPhysicsList * physicsList;
    G4VUserActionInitialization * userActionInitialization;
    G4UserWorkerInitialization * userWorkerInitialization;
    G4UserWorkerThreadInitialization * userWorkerThreadInitialization;
    G4UserRunAction * userRunAction;
    G4VUserPrimaryGeneratorAction * userPrimaryGeneratorAction;
    G4UserEventAction * userEventAction;
    G4UserStackingAction * userStackingAction;
    G4UserTrackingAction * userTrackingAction;
    G4UserSteppingAction * userSteppingAction;

  private:
    G4RunMessenger* runMessenger;

  protected:
    G4bool geometryInitialized;
    G4bool physicsInitialized;
    G4bool runAborted;
    G4bool initializedAtLeastOnce;
    G4bool geometryToBeOptimized;

    G4int runIDCounter;
    G4int verboseLevel;
    G4int printModulo;
    G4Timer * timer;
    G4DCtable* DCtable;

    G4Run* currentRun;
    G4Event* currentEvent;
    std::list<G4Event*>* previousEvents;
    G4int n_perviousEventsToBeStored;
    G4int numberOfEventToBeProcessed;

    G4bool storeRandomNumberStatus;
    G4int storeRandomNumberStatusToG4Event;
    G4String randomNumberStatusDir;
    G4String randomNumberStatusForThisRun;
    G4String randomNumberStatusForThisEvent;
    G4bool rngStatusEventsFlag;
    virtual void StoreRNGStatus(const G4String& filenamePrefix );
    
    G4VPhysicalVolume* currentWorld;

    G4int nParallelWorlds;

    G4String msgText;
    G4int n_select_msg;
    G4int numberOfEventProcessed;
    G4String selectMacro;
    G4bool fakeRun;

  public:
    virtual void rndmSaveThisRun();
    virtual void rndmSaveThisEvent();
    virtual void RestoreRandomNumberStatus(const G4String& fileN);

  public: // with description
    //The following set user-actions and user-initialization to the kernel
    //In MT mode, actions are shared among all threads, and should be set
    //in the master thread, while user-actions are thread-private and each      `
    //thread has private instances. Master thread does not have user-actions
    //except for the (optional) run-action.
    //User should instantiate the user-actions in the action-initialization
    //and use that class set method to set user-actions and not directly
    //the methods provided here.
    //Multiple Run,Event,Tracking, and Stepping actions are allowed, set
    //multiple instances and these will be appended to the current configuration
    //Multiple Stacking and PrimaryGeneration are not allowed
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

    //  These methods store respective user initialization and action classes.
    inline const G4VUserDetectorConstruction* GetUserDetectorConstruction() const
    { return userDetector; }
    inline const G4VUserPhysicsList* GetUserPhysicsList() const
    { return physicsList; }
    inline const G4VUserActionInitialization* GetUserActionInitialization() const
    { return userActionInitialization; }
    inline G4VUserActionInitialization* GetNonConstUserActionInitialization() const
    { return userActionInitialization; }
    inline const G4UserWorkerInitialization* GetUserWorkerInitialization() const
    { return userWorkerInitialization; }
    inline const G4UserWorkerThreadInitialization* GetUserWorkerThreadInitialization() const
    { return userWorkerThreadInitialization; }
    inline const G4UserRunAction* GetUserRunAction() const
    { return userRunAction; }
    inline const G4VUserPrimaryGeneratorAction* GetUserPrimaryGeneratorAction() const
    { return userPrimaryGeneratorAction; }
    inline const G4UserEventAction* GetUserEventAction() const
    { return userEventAction; }
    inline const G4UserStackingAction* GetUserStackingAction() const
    { return userStackingAction; }
    inline const G4UserTrackingAction* GetUserTrackingAction() const
    { return userTrackingAction; }
    inline const G4UserSteppingAction* GetUserSteppingAction() const
    { return userSteppingAction; }
    //  These methods returns respective user initialization and action classes.
     
    inline void SetNumberOfAdditionalWaitingStacks(G4int iAdd)
    { eventManager->SetNumberOfAdditionalWaitingStacks(iAdd); }
      //  Set the number of additional (optional) waiting stacks.
      // This method must be invoked at PreInit, Init or Idle states.
      // Once the user set the number of additional waiting stacks,
      // he/she can use the corresponding ENUM in G4ClassificationOfNewTrack.

    inline const G4String& GetVersionString() const
    { return kernel->GetVersionString(); }

    inline void SetPrimaryTransformer(G4PrimaryTransformer* pt)
    { kernel->SetPrimaryTransformer(pt); }

    inline void StoreRandomNumberStatusToG4Event(G4int vl)
      // if vl = 1 : status before primary particle generation is stored
      // if vl = 2 : status before event processing (after primary particle generation) is stored
      // if vl = 3 : both are stored
      // if vl = 0 : none is stored (default)
    { 
      storeRandomNumberStatusToG4Event = vl;
      eventManager->StoreRandomNumberStatusToG4Event(vl);
    }
    inline G4int GetFlagRandomNumberStatusToG4Event() const
    { return storeRandomNumberStatusToG4Event; }

  public:
    inline void SetRandomNumberStore(G4bool flag)
    { storeRandomNumberStatus = flag; }
    inline G4bool GetRandomNumberStore() const
    { return storeRandomNumberStatus; }
    inline void SetRandomNumberStoreDir(const G4String& dir)
    { 
      G4String dirStr = dir;
      if( dirStr(dirStr.length()-1) != '/' ) dirStr += "/";
#ifndef WIN32
      G4String shellCmd = "mkdir -p ";
#else
      std::replace(dirStr.begin(), dirStr.end(),'/','\\');
      G4String shellCmd = "if not exist " + dirStr + " mkdir ";
#endif
      shellCmd += dirStr;
      randomNumberStatusDir = dirStr;
      G4int sysret = system(shellCmd);
      if(sysret!=0)
      { 
        G4String errmsg = "\"" + shellCmd + "\" returns non-zero value. Directory creation failed.";
        G4Exception("GrRunManager::SetRandomNumberStoreDir","Run0071",JustWarning,errmsg);
        G4cerr << " return value = " << sysret << G4endl;
      }
    }
    inline const G4String& GetRandomNumberStoreDir() const
    { return randomNumberStatusDir; }
    inline const G4String& GetRandomNumberStatusForThisRun() const
    { return randomNumberStatusForThisRun; }
    inline const G4String& GetRandomNumberStatusForThisEvent() const
    {
      if(storeRandomNumberStatusToG4Event==0 || storeRandomNumberStatusToG4Event==2)
      { G4Exception("GrRunManager::SetRandomNumberStoreDir",
                    "Run0072",JustWarning,
                    "Random number status is not available for this event."); }
      return randomNumberStatusForThisEvent;
    }
    inline void SetRandomNumberStorePerEvent( G4bool flag )
    { rngStatusEventsFlag = flag; }
    inline G4bool GetRandomNumberStorePerEvent() const
    { return rngStatusEventsFlag; }
  public: // with description
    void GeometryHasBeenModified(G4bool prop=true);
    //  This method must be invoked (or equivalent UI command can be used)
    // in case the user changes his/her detector geometry after Initialize()
    // method has been invoked. Then, at the begining of the next BeamOn(),
    // all necessary re-voxelization will be made.
    //  The parameter "prop" has to be true if this C++ method is directly
    // invoked.

    void ReinitializeGeometry(G4bool destroyFirst=false, G4bool prop=true);
    //  This method must be invoked (or equivalent UI command can be used)
    // in case the user needs his/her detector construction has to be
    // re-invoked. Re-voxelization will be also done.
    //  If the first parameter "destroyFirst" is true, G4SolidStore,
    // G4LogicalVolumeStore and G4PhysicalVolumeStore are cleaned up, and
    // thus all solids, logical volumes and physical volumes previously defined
    // are deleted.
    //  The second parameter "prop" has to be true if this C++ method is directly
    // invoked.

    inline void PhysicsHasBeenModified()
    { kernel->PhysicsHasBeenModified(); }
    //  This method must be invoked (or equivalent UI command can be used)
    // in case the user changes his/her physics process(es), e.g. (in)activate 
    // some processes. Once this method is invoked, regardless of cuts are 
    // changed or not, BuildPhysicsTable() of PhysicsList is invoked for 
    // refreshing all physics tables.

    inline void CutOffHasBeenModified()
    {
      G4cerr << "CutOffHasBeenModified becomes obsolete." << G4endl;
      G4cerr << "It is safe to remove invoking this method." << G4endl;
    }  

  public: // with description
    void ReOptimizeMotherOf(G4VPhysicalVolume*);
    //  This method may be used if the orientation and/or size of this
    // particular physical volume has been modified while rest of the
    // geometries in the world has not been changed. This avoids the
    // full re-optimization of the entire geometry tree which is forced
    // if GeometryHasBeenModified() method is invoked.

    void ReOptimize(G4LogicalVolume*);
    //  Same as above, but the mother logical volume is specified.

  public:
    inline void SetVerboseLevel(G4int vl)
    { verboseLevel = vl; 
      kernel->SetVerboseLevel(vl); }
    inline G4int GetVerboseLevel() const
    { return verboseLevel; }
    inline G4int GetPrintProgress()
    { return printModulo; }
    inline void SetPrintProgress(G4int i)
    { printModulo = i; }

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
    { return geometryToBeOptimized; }

  public: // with description
    inline void SetNumberOfEventsToBeStored(G4int val)
    { n_perviousEventsToBeStored = val; }
    //  Sets the number of events to be kept after processing. That is, "val" previous
    // events can be used with the most recent event for digitizing pileup. "val"+1
    // previous event is deleted.
    //  This method must be invoked before starting the event loop.
    inline const G4Run* GetCurrentRun() const
    { return currentRun; }
    inline G4Run* GetNonConstCurrentRun() const
    { return currentRun; }
    //  Returns the pointer to the current run. This method is available for Geant4
    // states of GeomClosed and EventProc.
    inline const G4Event* GetCurrentEvent() const
    { return currentEvent; }
    //  Returns the pointer to the current event. This method is available for EventProc
    // state.
    inline const G4Event* GetPreviousEvent(G4int i) const
    {
      if(i>=1 && i<=n_perviousEventsToBeStored)
      {
        std::list<G4Event*>::iterator itr = previousEvents->begin();
        for(G4int j=1;j<i;j++) { itr++; }
        return *itr;
      }
      return 0;
    }
    //  Returns the pointer to the "i" previous event. This method is availavle for
    // EventProc state. In case the event loop has not yet to reach to the requested
    // event, null will be returned. To use this method, SetNumberOfEventsToBeStored()
    // method mentioned above must be invoked previously to the event loop.
    inline void SetRunIDCounter(G4int i)
    { runIDCounter = i; }
    //  Set the run number counter. Initially, the counter is initialized to zero and
    // incremented by one for every BeamOn().

  public:
    inline G4int GetNumberOfParallelWorld() const
    { return nParallelWorlds; }
    inline void SetNumberOfEventsToBeProcessed(G4int val)
    { numberOfEventToBeProcessed = val; }
    inline G4int GetNumberOfEventsToBeProcessed() const
    { return numberOfEventToBeProcessed; }
    inline G4int GetNumberOfSelectEvents() const
    { return n_select_msg; }
    inline G4String GetSelectMacro() const
    { return selectMacro; }
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }

  public:
    inline RMType GetRunManagerType() const
    { return runManagerType; }

  protected:
    RMType runManagerType;

  public:
    virtual void ConstructScoringWorlds();
  protected:
    void UpdateScoring();
    virtual void DeleteUserInitializations();
    //Called by destructor to delete user detector. Note: the userdetector is shared by threads
    //Thus this should be re-implemented to empty in derived classes that implement the worker model
  private:
    //disable assignment and copy constructors
    G4RunManager(const G4RunManager&) {}
    G4RunManager& operator=(const G4RunManager&) { return *this; }

  protected:
    // This boolean flag has to be shared by all G4RunManager objects
    static G4bool fGeometryHasBeenDestroyed;
  public:
    static G4bool IfGeometryHasBeenDestroyed();
    //This is used only by workers thread to reset RNG engines from files
    //that are event specific. Not implemented for sequential since run seed
    //defines event seeds
    virtual void RestoreRndmEachEvent(G4bool) { /*No effect in SEQ */ }
};

#endif

