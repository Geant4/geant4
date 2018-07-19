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
// $Id: G4EventManager.hh 94950 2016-01-07 11:53:14Z gcosmo $
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
class G4StateManager;
#include "globals.hh"
class G4VUserEventInformation;

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
      static G4ThreadLocal G4EventManager* fpEventManager;

  public:
      G4EventManager();
      ~G4EventManager();

  private:
      G4EventManager(const G4EventManager &right) = delete;
      G4EventManager& operator=(const G4EventManager& right) = delete;

  public: // with description
      void ProcessOneEvent(G4Event* anEvent);
      //  This method is the main entry to this class for simulating an event.

      void ProcessOneEvent(G4TrackVector* trackVector,G4Event* anEvent=0);
      //  This is an alternative entry for large HEP experiments which create G4Track
      // objects by themselves directly without using G4VPrimaryGenerator or user
      // primary generator action. Dummy G4Event object will be created if "anEvent" is null
      // for internal use, but this dummy object will be deleted at the end of this
      // method and will never be available for the use after the processing.
      // Note that in this case of null G4Event pointer no output of the simulated event 
      // is returned by this method, but the user must implement some mechanism
      // of storing output by his/herself, e.g. in his/her UserEventAction and/or
      // sensitive detectors.
      //  If valid G4Event object is given, this object will not be deleted with
      // this method and output objects such as hits collections and trajectories
      // will be associated to this event object. If this event object has valid
      // primary vertices/particles, they will be added to the given trackvector input.

  private:
      void DoProcessing(G4Event* anEvent);
      void StackTracks(G4TrackVector *trackVector, G4bool IDhasAlreadySet=false);
  
      G4Event* currentEvent;

      G4StackManager *trackContainer;
      G4TrackingManager *trackManager;
      G4TrajectoryContainer *trajectoryContainer;
      G4int trackIDCounter;
      G4int verboseLevel;
      G4SDManager* sdManager;
      G4PrimaryTransformer* transformer;
      G4bool tracking;
      G4bool abortRequested;

      G4EvManMessenger* theMessenger;

      G4UserEventAction*    userEventAction;
      G4UserStackingAction* userStackingAction;
      G4UserTrackingAction* userTrackingAction;
      G4UserSteppingAction* userSteppingAction;

      G4int storetRandomNumberStatusToG4Event;
      G4String randomNumberStatusToG4Event;

      G4StateManager* stateManager;

  public: // with description
      inline const G4Event* GetConstCurrentEvent()
      { return currentEvent; }
      inline G4Event* GetNonconstCurrentEvent()
      { return currentEvent; }
      //  These methods returns the pointers of const G4Event*
      // and G4Event*, respectively. Null will be returned when
      // an event is not processing.

  public: // with description
      void AbortCurrentEvent();
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
      void SetNumberOfAdditionalWaitingStacks(G4int iAdd)
      { trackContainer->SetNumberOfAdditionalWaitingStacks(iAdd); }

      void KeepTheCurrentEvent();
      // If the current event exists, it is kept undeleted until the end of the current run

      inline G4StackManager* GetStackManager() const
      { return trackContainer; }
      inline G4TrackingManager* GetTrackingManager() const
      { return trackManager; }

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

      void SetUserInformation(G4VUserEventInformation* anInfo);
      G4VUserEventInformation* GetUserInformation();
      // Set and get method of G4VUserEventInformation object associating with
      // the current event. Both methods are valid only for G4State_EventProc
      // application state.

      inline G4PrimaryTransformer* GetPrimaryTransformer() const
      { return transformer; }
      inline void SetPrimaryTransformer(G4PrimaryTransformer* tf)
      { transformer = tf; }
      inline void StoreRandomNumberStatusToG4Event(G4int vl)
      { storetRandomNumberStatusToG4Event = vl; }
};



#endif

