// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RunManager.hh,v 1.1 1999-01-07 16:14:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a class for run controle in GEANT4
// 
//      User must provide his own classes derived from the following
//      three abstract classes and register them to the RunManager. 
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
//	History
//        first version                              by Makoto Asai 
//        modified                       11 Mar. 1998 by H.Kurashige
// ------------------------------------------------------------

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
#include <rw/tpordvec.h>

class G4RunManager
{
  public:
    static G4RunManager* GetRunManager();

  private:
    static G4RunManager* fRunManager;

  public:
    G4RunManager();
    virtual ~G4RunManager();

  public:
    virtual void BeamOn(G4int n_event,const char* macroFile=NULL,G4int n_select=-1);
    virtual void Initialize();
    virtual void DefineWorldVolume(G4VPhysicalVolume * worldVol);
    virtual void AbortRun();

  protected:
    // methods invoked from Initialize()
    virtual void InitializeGeometry();
    virtual void InitializePhysics();
    virtual void InitializeCutOff();

    // methods invoked from BeamOn()
    virtual G4bool ConfirmBeamOnCondition();
    virtual void RunInitialization();
    virtual void DoEventLoop(G4int n_event,const char* macroFile=NULL,G4int n_select=-1);
    virtual void RunTermination();

    // methods invoked from DoEventLoop()
    virtual G4Event* GenerateEvent(G4int i_event);
    virtual void AnalyzeEvent(G4Event* anEvent);

  protected:
    void StackPreviousEvent(G4Event* anEvent);

  protected:
    G4EventManager * eventManager;
    G4VUserDetectorConstruction * userDetector;
    G4VUserPhysicsList * physicsList;
    G4UserRunAction * userRunAction;
    G4VUserPrimaryGeneratorAction * userPrimaryGeneratorAction;

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

    G4int verboseLevel;
    G4Timer * timer;
    G4DCtable* DCtable;

    G4Run* currentRun;
    G4Event* currentEvent;
    RWTPtrOrderedVector<G4Event>* previousEvents;
    G4int n_perviousEventsToBeStored;

  public:
    inline void SetUserInitialization(G4VUserDetectorConstruction* userInit)
    { userDetector = userInit; }
    inline void SetUserInitialization(G4VUserPhysicsList* userInit)
    { physicsList = userInit; }
    inline void SetUserAction(G4UserRunAction* userAction)
    { userRunAction = userAction; }
    inline void SetUserAction(G4VUserPrimaryGeneratorAction* userAction)
    { userPrimaryGeneratorAction = userAction; }
    inline void SetUserAction(G4UserEventAction* userAction)
    { eventManager->SetUserAction(userAction); }
    inline void SetUserAction(G4UserStackingAction* userAction)
    { eventManager->SetUserAction(userAction); }
    inline void SetUserAction(G4UserTrackingAction* userAction)
    { eventManager->SetUserAction(userAction); }
    inline void SetUserAction(G4UserSteppingAction* userAction)
    { eventManager->SetUserAction(userAction); }

    inline void GeometryHasBeenModified()
    { geometryNeedsToBeClosed = true; }
    inline void CutOffHasBeenModified()
    { cutoffInitialized = false; }

    inline void SetPauseAtBeginOfEvent(G4bool vl)
    { pauseAtBeginOfEvent = vl; }
    inline void SetPauseAtEndOfEvent(G4bool vl)
    { pauseAtEndOfEvent = vl; }

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

    inline void SetNumberOfEventsToBeStored(G4int val)
    { n_perviousEventsToBeStored = val; }
    inline const G4Run* GetCurrentRun() const
    { return currentRun; }
    inline const G4Event* GetCurrentEvent() const
    { return currentEvent; }
    inline const G4Event* GetPreviousEvent(G4int i) const
    {
      if(i>=1 && i<=n_perviousEventsToBeStored)
      { return (*previousEvents)[i-1]; }
      return NULL;
    }

    inline void SetDCtable(G4DCtable* DCtbl)
    { DCtable = DCtbl; }
};

#endif

