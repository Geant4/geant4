//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EventManager.hh,v 1.12 2003/09/09 20:09:17 asaim Exp $
// GEANT4 tag $Name: geant4-06-00 $
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
class G4VUserEventInformation;

#ifndef WIN32         // Temporarly disabled on Windows, until CLHEP
                      // will support the HepMC module
#include "CLHEP/HepMC/GenEvent.h"
#endif

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
      //  This method is the main entry to this class for simulating an event.

#ifndef WIN32         // Temporarly disabled on Windows, until CLHEP
                      // will support the HepMC module
      void ProcessOneEvent(const HepMC::GenEvent* hepmcevt,G4Event* anEvent=0);
      //  This is an alternative entry for large HEP experiments which use
      // HepMC event class. Dummy G4Event object will be created if "anEvent" is null
      // for internal use, but this dummy object will be deleted at the end of this
      // method and will never be available for the use after the processing.
      // Note that in this case of null G4Event pointer no output of the simulated event
      // is returned by this method, but the user must implement some mechanism
      // of storing output by his/herself, e.g. in his/her UserEventAction and/or
      // sensitive detectors.
      //  If valid G4Event object is given, this object will not be deleted with
      // this method and output objects such as hits collections and trajectories
      // will be associated to this event object. If this event object has valid
      // primary vertices/particles, they will be added to the given HepMC event input.
#endif

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
      void SetNumberOfAdditionalWaitingStacks(G4int iAdd)
      { trackContainer->SetNumberOfAdditionalWaitingStacks(iAdd); }

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

};



#endif

