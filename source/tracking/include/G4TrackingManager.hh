// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackingManager.hh,v 1.1 1999-01-07 16:14:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4TrackingManager.hh
//
// Description:
//   This class represents the manager who steers the tracking
//   of a particle given from the event manager.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

class G4TrackingManager;

#ifndef G4TrackingManager_h
#define G4TrackingManager_h 1

#include "globals.hh"                  // Include from 'global'
/////#include "G4Hit.hh"               // Include from 'Hit/dig'
#include "G4SteppingManager.hh"        // Include from 'tracking'
#include "G4Track.hh"                  // Include from 'tracking'
#include "G4TrackVector.hh"            // Include from 'tracking'
#include "G4TrackStatus.hh"            // Include from 'tracking'
#include "G4StepStatus.hh"             // Include from 'tracking'
#include "G4UserTrackingAction.hh"     // Include from 'tracking'
#include "G4UserSteppingAction.hh"     // Include from 'tracking'
#include "G4Trajectory.hh"             // Include from 'tracking'


////////////////////////
class G4TrackingManager 
////////////////////////
{

//--------
   public:
//--------

// Constructor/Destructor

   G4TrackingManager();
      // TrackingManger should be dynamic persistent, therefore you
      // need to invoke new() when you call this constructor.
      // "G4SteppingManger' and "G4UserTrackingAction" will be 
      // constructed in this constructor. "This" pointer will be
      // passed to "G4UserTrackingAction". 

   ~G4TrackingManager();

// Get/Set functions

   G4Track* GetTrack() const;

   G4bool GetStoreTrajectory() const;
   void SetStoreTrajectory(G4bool value);

   G4SteppingManager* GetSteppingManager() const;

   G4UserTrackingAction* GetUserTrackingAction() const;

   G4Trajectory* GimmeTrajectory() const;
    
   G4TrackVector* GimmeSecondaries() const;

   void SetNavigator(G4Navigator* apValue);

   void SetUserAction(G4UserTrackingAction* apAction);
   void SetUserAction(G4UserSteppingAction* apAction);

   void SetVerboseLevel(G4int vLevel);
   G4int GetVerboseLevel() const;

// Other member functions

   void ProcessOneTrack(G4Track* apValueG4Track);
      // Invoking this function, a G4Track given by the argument
      // will be tracked.  

   void EventAborted();
      // Invoking this function, the current tracking will be
      // aborted immediately. The tracking will return the 
      // G4TrackStatus in 'fUserKillTrackAndSecondaries'.
      // By this the EventManager deletes the current track and all 
      // its accoicated csecondaries.

//---------
   private:
//---------

// Member functions

   void Verbose(G4String select);
      // Output messages according to the current value of the 
      // Verbose level.

// Member data

   G4Track* fpTrack;
   G4SteppingManager* fpSteppingManager;
   G4UserTrackingAction* fpUserTrackingAction;
   G4Trajectory* fpTrajectory;
   G4bool StoreTrajectory;
   G4int verboseLevel;

};


//*******************************************************************
//
//  Inline function 
//
//*******************************************************************

   inline G4Track* G4TrackingManager::GetTrack() const { 
     return fpTrack;
   }

   inline G4bool G4TrackingManager::GetStoreTrajectory() const { 
     return StoreTrajectory;
   }

   inline void G4TrackingManager::SetStoreTrajectory(G4bool value){ 
     StoreTrajectory = value;
   }

   inline G4SteppingManager* G4TrackingManager::GetSteppingManager() const { 
     return fpSteppingManager; 
   }

   inline G4UserTrackingAction* G4TrackingManager::GetUserTrackingAction() const  { 
     return fpUserTrackingAction; 
   }

   inline G4Trajectory* G4TrackingManager::GimmeTrajectory() const { 
     return fpTrajectory ; 
   }
    
   inline G4TrackVector* G4TrackingManager::GimmeSecondaries() const { 
     return fpSteppingManager->GetSecondary(); 
   }

   inline void G4TrackingManager::SetUserAction(G4UserTrackingAction* apAction){
     if (fpUserTrackingAction) delete fpUserTrackingAction;
     fpUserTrackingAction = apAction;
     apAction->SetTrackingManagerPointer(this);
   }

   inline void G4TrackingManager::SetUserAction(G4UserSteppingAction* apAction){
     fpSteppingManager->SetUserAction(apAction);
     apAction->SetSteppingManagerPointer(fpSteppingManager);  
   }

   inline void G4TrackingManager::SetVerboseLevel(G4int vLevel){ 
     verboseLevel = vLevel; 
     fpSteppingManager -> SetVerboseLevel( vLevel );
   }


   inline G4int G4TrackingManager::GetVerboseLevel() const { 
     return verboseLevel; 
   }


#endif
