// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EventManager.hh,v 1.5 2000-01-12 01:29:50 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//


#ifndef G4EventManager_h
#define G4EventManager_h 1

#include "evmandefs.hh"
#include "G4StackManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4PrimaryTransformer.hh"
class G4Event;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4EvManMessenger;
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4VTrajectory.hh"
#include "G4TrackStatus.hh"
class G4SDManager;
#include "globals.hh"

// class description:
//
//	G4EventManager controls an event. This class must be a singleton
//      and should be constructed by G4RunManager.
//

class G4EventManager 
{
  public: // with description
      static G4EventManager* GetEventManager();
      //  This method returns the singleton pointer of G4EventManager.

  private:
      static G4EventManager* fpEventManager;

  public:
      G4EventManager();
      ~G4EventManager();

  private:
      G4EventManager(const G4EventManager &right);
      G4EventManager& operator=(const G4EventManager& right);

  public: // with description
      void ProcessOneEvent(G4Event* anEvent);
      //  This method it the main entry to this class for simulating an event.
      // This method must be exclusively invoked by G4RunManager.

  private:
      void StackTracks(G4TrackVector *trackVector);
  
      G4Event* currentEvent;

      G4StackManager *trackContainer;
      G4TrackingManager *trackManager;
      G4TrajectoryContainer *trajectoryContainer;
      G4int trackIDCounter;
      int verboseLevel;
      G4SDManager* sdManager;
      G4PrimaryTransformer* transformer;
      G4bool tracking;

      G4EvManMessenger* theMessenger;

      G4UserEventAction*    userEventAction;
      G4UserStackingAction* userStackingAction;
      G4UserTrackingAction* userTrackingAction;
      G4UserSteppingAction* userSteppingAction;

  public: // with description
      inline const G4Event* GetConstCurrentEvent()
      { return currentEvent; }
      inline G4Event* GetNonconstCurrentEvent()
      { return currentEvent; }
      //  These methods returns the pointers of const G4Event*
      // and G4Event*, respectively. Null will be returned when
      // an event is not processing.

  public: // with description
      inline void AbortCurrentEvent()
      { 
        trackContainer->clear();
        if(tracking) trackManager->EventAborted();
      }
      //  This method aborts the processing of the current event. All stacked
      // tracks are deleted. The contents of G4Event object is not completed,
      // but trajectories, hits, and/or digits which are created before the
      // moment of abortion can be used.

  public: // with description
      void SetUserAction(G4UserEventAction* userAction);
      void SetUserAction(G4UserStackingAction* userAction);
      void SetUserAction(G4UserTrackingAction* userAction);
      void SetUserAction(G4UserSteppingAction* userAction);
      inline G4UserEventAction* GetUserEventAction()
      { return userEventAction; }
      inline G4UserStackingAction* GetUserStackingAction()
      { return userStackingAction; }
      inline G4UserTrackingAction* GetUserTrackingAction()
      { return userTrackingAction; }
      inline G4UserSteppingAction* GetUserSteppingAction()
      { return userSteppingAction; }
      // Set and get methods for user action classes. User action classes
      // which should belong to the other managers will be sent to the 
      // corresponding managers.

  public: // with description
      inline G4int GetVerboseLevel()
      { return verboseLevel; }
      inline void SetVerboseLevel( G4int value )
      {
        verboseLevel = value;
        trackContainer->SetVerboseLevel( value );
        transformer->SetVerboseLevel( value );
      }
      // Set and get method of the verbose level

};



#endif

