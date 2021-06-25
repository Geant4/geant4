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
// G4TrackingManager class implementation
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------

#include "G4TrackingManager.hh"
#include "G4Trajectory.hh"
#include "G4SmoothTrajectory.hh"
#include "G4RichTrajectory.hh"
#include "G4ios.hh"
#include "G4Profiler.hh"
#include "G4TiMemory.hh"

//////////////////////////////////////
G4TrackingManager::G4TrackingManager()
//////////////////////////////////////
{
  messenger = new G4TrackingMessenger(this);
  fpSteppingManager = new G4SteppingManager();
}

///////////////////////////////////////
G4TrackingManager::~G4TrackingManager()
///////////////////////////////////////
{
  delete messenger;
  delete fpSteppingManager;
  delete fpUserTrackingAction;
}

////////////////////////////////////////////////////////////////
void G4TrackingManager::ProcessOneTrack(G4Track* apValueG4Track)
////////////////////////////////////////////////////////////////
{
  // Receiving a G4Track from the EventManager, this funciton has the
  // responsibility to trace the track till it stops

  fpTrack = apValueG4Track;
  EventIsAborted = false;

  // Clear secondary particle vector
  //
  for(std:: size_t itr=0; itr<GimmeSecondaries()->size(); ++itr)
  { 
    delete (*GimmeSecondaries())[itr];
  }
  GimmeSecondaries()->clear();  
   
  if(verboseLevel>0 && (G4VSteppingVerbose::GetSilent()!=1) ) TrackBanner();
  
  // Give SteppingManger the pointer to the track which will be tracked
  //
  fpSteppingManager->SetInitialStep(fpTrack);

  // Pre tracking user intervention process

  fpTrajectory = nullptr;
  if( fpUserTrackingAction != nullptr )
  {
    fpUserTrackingAction->PreUserTrackingAction(fpTrack);
  }

  // we need this to scope the G4Track::ProfilerConfig b/t
  // the PreUserTrackingAction and PostUserTrackingAction
  {
#if defined(GEANT4_USE_TIMEMORY)
    ProfilerConfig profiler{ fpTrack };
#endif

#ifdef G4_STORE_TRAJECTORY
    // Construct a trajectory if it is requested
    //
    if(StoreTrajectory && (fpTrajectory == nullptr))
    {
      // default trajectory concrete class object
      switch(StoreTrajectory)
      {
        default:
        case 1:
          fpTrajectory = new G4Trajectory(fpTrack);
          break;
        case 2:
          fpTrajectory = new G4SmoothTrajectory(fpTrack);
          break;
        case 3:
          fpTrajectory = new G4RichTrajectory(fpTrack);
          break;
        case 4:
          fpTrajectory = new G4RichTrajectory(fpTrack);
          break;
      }
    }
#endif

    // Give SteppingManger the maxmimum number of processes
    fpSteppingManager->GetProcessNumber();

    // Give track the pointer to the Step
    fpTrack->SetStep(fpSteppingManager->GetStep());

    // Inform beginning of tracking to physics processes
    fpTrack->GetDefinition()->GetProcessManager()->StartTracking(fpTrack);

    // Track the particle Step-by-Step while it is alive
    //
    while((fpTrack->GetTrackStatus() == fAlive) ||
          (fpTrack->GetTrackStatus() == fStopButAlive))
    {
      fpTrack->IncrementCurrentStepNumber();
      fpSteppingManager->Stepping();
#ifdef G4_STORE_TRAJECTORY
      if(StoreTrajectory)
      {
        fpTrajectory->AppendStep(fpSteppingManager->GetStep());
      }
#endif
      if(EventIsAborted)
      {
        fpTrack->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
    // Inform end of tracking to physics processes
    fpTrack->GetDefinition()->GetProcessManager()->EndTracking();
  }

  // Post tracking user intervention process.
  if( fpUserTrackingAction != nullptr )
  {
    fpUserTrackingAction->PostUserTrackingAction(fpTrack);
  }

  // Destruct the trajectory if it was created
#ifdef G4VERBOSE
  if(StoreTrajectory&&verboseLevel>10)
  {
    fpTrajectory->ShowTrajectory();
  }
#endif
  if( (!StoreTrajectory) && (fpTrajectory != nullptr) )
  {
    delete fpTrajectory;
    fpTrajectory = nullptr;
  }
}

//////////////////////////////////////
void G4TrackingManager::SetTrajectory(G4VTrajectory* aTrajectory)
//////////////////////////////////////
{
#ifndef G4_STORE_TRAJECTORY
  G4Exception("G4TrackingManager::SetTrajectory()",
              "Tracking0015", FatalException,
              "Invoked without G4_STORE_TRAJECTORY option set!");
#endif
  fpTrajectory = aTrajectory;
}

//////////////////////////////////////
void G4TrackingManager::EventAborted()
//////////////////////////////////////
{
  fpTrack->SetTrackStatus( fKillTrackAndSecondaries );
  EventIsAborted = true;
}

//////////////////////////////////////
void G4TrackingManager::TrackBanner()
//////////////////////////////////////
{
  G4cout << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << "* G4Track Information: "
         << "  Particle = " << fpTrack->GetDefinition()->GetParticleName()
         << ","
         << "   Track ID = " << fpTrack->GetTrackID()
         << ","
         << "   Parent ID = " << fpTrack->GetParentID()
         << G4endl;
  G4cout << "*******************************************************"
         << "**************************************************"
         << G4endl;
  G4cout << G4endl;
}
