// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RunManager.hh,v 1.7 1999-11-11 15:37:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// class description:
//
//      This is a class for run control in GEANT4
// 
//     User must provide his own classes derived from the following
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
//     G4RunManager is the only manager class in Geant4 kernel which 
//     the user MUST construct an object by him/herself in the main(). 
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
class G4UserRunAction;
class G4VUserPrimaryGeneratorAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;

class G4VPhysicalVolume;
class G4Timer;
class G4RunMessenger;
class G4DCtable;
class G4Run;
#include "G4Event.hh"

#include "G4EventManager.hh"
#include "globals.hh"
#include "g4rw/tpordvec.h"

class G4RunManager
{
  public: // with description
    static G4RunManager* GetRunManager();
    //  Static method which returns the singleton pointer of G4RunManager or
    // its derived class.

  private:
    static G4RunManager* fRunManager;

  public: // with description
    G4RunManager();
    virtual ~G4RunManager();
    //  The constructor and the destructor. The user must construct this class
    // object at the beginning of his/her main() and must delete it at the 
    // bottom of the main().

  public: // with description
    virtual void BeamOn(G4int n_event,const char* macroFile=NULL,G4int n_select=-1);
    //  This method starts an event loof of "n_event" events. The condition of Geant4
    // is examined before starting the event loop. This method must be invoked at
    // Idle state. The state will be changed to GeomClosed during the event loop and
    // will go back to Idle when the loop is over or aborted.
    //  In case a string "macroFile" which represents the name of a macro file is given,
    // this macro file will be executed AT THE END of each event processing. In case
    // "n_select" is greater than zero, at the ond of first "n_select" events the macro
    // file is executed.
    virtual void Initialize();
    //  This method invokes all the necessary initialization procedures for an event
    // loop. This method must be invoked at the Geant4 state of PreInit or Idle. The
    // state will be changed to Init during the initialization procedures and then
    // changed to Idle.
    //  This method invokes three protected methods, InitializeGeometry(), 
    // InitializePhysics(), and InitializeCutOff().
    //  After some event loops, the user can invoke this method once again. It is
    // required if the user changes geometry, physics process, and/or cut off value.
    // If the user forget the second invokation, G4RunManager will invoke BeamOn()
    // method will invoke this method. (Note that this feature is not valid for the
    // first initialization.)
    virtual void DefineWorldVolume(G4VPhysicalVolume * worldVol);
    //  Usually, this method is invoked from InitializeGeometry() protected method
    // of this class. But, in case all of geometry has already created and kept in
    // the ODBMS, the pointer to the world physical volume can be set by this method.
    virtual void AbortRun();
    //  This method safely aborts the current event loop even if an event is in progress.
    // This method is available for Geant4 states of GeomClosed and EventProc. The state
    // will be changed to Idle, so that another event loop can be done.

  protected: // with description

    virtual void InitializeGeometry();
    virtual void InitializePhysics();
    virtual void InitializeCutOff();
    //  These three protected methods are invoked from Initialize() method for the 
    // initializations of geometry, physics processes, and cut off. The user's concrete
    // G4VUserDetectorConstruction class will be accessed from InitializeGeometry() and
    // G4VUserPhysicsList class will be accessed from other two methods.

    virtual G4bool ConfirmBeamOnCondition();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event,const char* macroFile=NULL,G4int n_select=-1);
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

  protected:
    void StackPreviousEvent(G4Event* anEvent);

  protected:
    G4EventManager * eventManager;
    G4VUserDetectorConstruction * userDetector;
    G4VUserPhysicsList * physicsList;
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
    G4bool cutoffInitialized;
    G4bool geometryNeedsToBeClosed;
    G4bool runAborted;
    G4bool initializedAtLeastOnce;
    G4bool geometryToBeOptimized;

    G4bool pauseAtBeginOfEvent;
    G4bool pauseAtEndOfEvent;

    G4int runIDCounter;
    G4int verboseLevel;
    G4Timer * timer;
    G4DCtable* DCtable;

    G4Run* currentRun;
    G4Event* currentEvent;
    G4RWTPtrOrderedVector<G4Event>* previousEvents;
    G4int n_perviousEventsToBeStored;

    G4int storeRandomNumberStatus;

  public:
    virtual void StoreRandomNumberStatus(G4int eventID=-1);
    virtual void RestoreRandomNumberStatus(G4String fileN);

  public: // with description
    inline void SetUserInitialization(G4VUserDetectorConstruction* userInit)
    { userDetector = userInit; }
    inline void SetUserInitialization(G4VUserPhysicsList* userInit)
    { physicsList = userInit; }
    inline void SetUserAction(G4UserRunAction* userAction)
    { userRunAction = userAction; }
    inline void SetUserAction(G4VUserPrimaryGeneratorAction* userAction)
    { userPrimaryGeneratorAction = userAction; }
    inline void SetUserAction(G4UserEventAction* userAction)
    { 
      eventManager->SetUserAction(userAction); 
      userEventAction = userAction;
    }
    inline void SetUserAction(G4UserStackingAction* userAction)
    { 
      eventManager->SetUserAction(userAction); 
      userStackingAction = userAction;
    }
    inline void SetUserAction(G4UserTrackingAction* userAction)
    { 
      eventManager->SetUserAction(userAction); 
      userTrackingAction = userAction;
    }
    inline void SetUserAction(G4UserSteppingAction* userAction)
    { 
      eventManager->SetUserAction(userAction); 
      userSteppingAction = userAction;
    }
    //  These methods store respective user initialization and action classes.
    inline const G4VUserDetectorConstruction* GetUserDetectorConstruction() const
    { return userDetector; }
    inline const G4VUserPhysicsList* GetUserPhysicsList() const
    { return physicsList; }
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

  public:
    inline void SetRandomNumberStore(G4int i)
    { storeRandomNumberStatus = i; }
    inline G4int GetRandomNumberStore() const
    { return storeRandomNumberStatus; }

  public: // with description
    inline void GeometryHasBeenModified()
    { geometryNeedsToBeClosed = true; }
    inline void CutOffHasBeenModified()
    { cutoffInitialized = false; }
    //  These two methods must be invoked (or equivalent UI commands can be used)
    // in case the user changes his/her detector geometry or cut off value(s) after
    // Initialize() metho has been invoked. Then, at the begining of the next BeamOn(),
    // all necessary re-initialization will be done.
    inline void SetPauseAtBeginOfEvent(G4bool vl)
    { pauseAtBeginOfEvent = vl; }
    inline void SetPauseAtEndOfEvent(G4bool vl)
    { pauseAtEndOfEvent = vl; }
    //  If the boolean flags are true, Pause() method of G4StateManager is invoked
    // at the very begining (before generating a G4Event object) or at the end of
    // each event. So that, in case a (G)UI session is defined, the user can interact.

  public:
    inline void SetVerboseLevel(G4int vl)
    { verboseLevel = vl; }
    inline G4int GetVerboseLevel() const
    { return verboseLevel; }

    inline void SetGeometryToBeOptimized(G4bool vl)
    { 
      if(geometryToBeOptimized != vl)
      {
        geometryToBeOptimized = vl;
        geometryNeedsToBeClosed = true;
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
    //  Returns the pointer to the current run. This method is available for Geant4
    // states of GeomClosed and EventProc.
    inline const G4Event* GetCurrentEvent() const
    { return currentEvent; }
    //  Returns the pointer to the current event. This method is available for EventProc
    // state.
    inline const G4Event* GetPreviousEvent(G4int i) const
    {
      if(i>=1 && i<=n_perviousEventsToBeStored)
      { return (*previousEvents)[i-1]; }
      return NULL;
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
    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
};

#endif

