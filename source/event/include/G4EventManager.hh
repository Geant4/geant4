// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4EventManager.hh,v 1.2 1999-04-15 08:41:50 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 10/Dec/96 M.Asai
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

//## Class: G4EventManager
//	G4EventManager controls an event. This class has only
//	one public method, ProcessOneEvent(), which will be
//	called by GEANT4.

class G4EventManager 
{

  public:
      G4EventManager();
      ~G4EventManager();

  private:
      G4EventManager(const G4EventManager &right);
      G4EventManager& operator=(const G4EventManager& right);

  public:
      void ProcessOneEvent(G4Event* anEvent);

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
      G4UserEventAction* userEventAction;
      G4bool tracking;

      G4EvManMessenger* theMessenger;

  public:
      inline const G4Event* GetConstCurrentEvent()
      { return currentEvent; };
      inline G4Event* GetNonconstCurrentEvent()
      { return currentEvent; };
      inline void AbortCurrentEvent()
      { 
        trackContainer->clear();
        if(tracking) trackManager->EventAborted();
      };
      void SetUserAction(G4UserEventAction* userAction);
      inline void SetUserAction(G4UserStackingAction* userAction)
      { trackContainer->SetUserStackingAction(userAction); };
      inline void SetUserAction(G4UserTrackingAction* userAction)
      { trackManager->SetUserAction(userAction); };
      inline void SetUserAction(G4UserSteppingAction* userAction)
      { trackManager->SetUserAction(userAction); };
      inline G4int GetVerboseLevel()
      { return verboseLevel; };
      inline void SetVerboseLevel( G4int value )
      {
        verboseLevel = value;
        trackContainer->SetVerboseLevel( value );
        transformer->SetVerboseLevel( value );
      };

};



#endif

