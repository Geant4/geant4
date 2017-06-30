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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
//
// $Id$
//
/// \file ITTrackingInteractivity.hh
/// \brief Implementation of the ITTrackingInteractivity class

#include "ITTrackingInteractivity.hh"
#include "G4TrackingInformation.hh"
#include "G4VTrajectory.hh"
#include "G4Trajectory.hh"
#include "G4SmoothTrajectory.hh"
#include "G4RichTrajectory.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4IT.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4VSteppingVerbose.hh"
#include "G4VisManager.hh"

class G4Trajectory_Lock
{
  friend class ITTrackingInteractivity;

  G4Trajectory_Lock(): fpTrajectory(0)
  {;}

  ~G4Trajectory_Lock()
  {;}

  G4VTrajectory* fpTrajectory;
};

ITTrackingInteractivity::ITTrackingInteractivity()
{
  fStoreTrajectory = 0;
  fVerboseLevel = 0;

  fpUserTrackingAction  = 0;
  fpUserSteppingAction = 0;

  ////////////////////////////
  // In case you want to use same tracking/stepping action
  // for normal and IT stepping
  /*
    fpUserTrackingAction =
      trackingManager->GetUserTrackingAction();
    fpUserSteppingAction =
      G4EventManager::GetEventManager()->GetUserSteppingAction();
   */
  ////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ITTrackingInteractivity::~ITTrackingInteractivity()
{
  G4EventManager* eventManager =
      G4EventManager::GetEventManager();

  if (eventManager)
  {
  G4UserTrackingAction* std_trackAct =
      eventManager->GetUserTrackingAction();
  if(fpUserTrackingAction != std_trackAct && fpUserTrackingAction)
    delete fpUserTrackingAction;

  G4UserSteppingAction* std_stepAct =
      eventManager->GetUserSteppingAction();
  if(fpUserSteppingAction != std_stepAct && fpUserSteppingAction)
    delete fpUserSteppingAction;
  }
  else
  {
    if(fpUserSteppingAction)
    {
      delete fpUserSteppingAction;
    }
    
    if(fpUserTrackingAction)
    {
      delete fpUserTrackingAction;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::Initialize()
{
  G4TrackingManager* trackingManager = 
    G4EventManager::GetEventManager()->GetTrackingManager();
  fStoreTrajectory = trackingManager->GetStoreTrajectory();
  fVerboseLevel = trackingManager->GetVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::StartTracking(G4Track* track)
{
#ifdef G4VERBOSE
  if(fVerboseLevel)
  {
    TrackBanner(track, "G4ITTrackingManager::StartTracking : ");
  }
#endif

  if(fVerboseLevel>0 && (G4VSteppingVerbose::GetSilent()!=1) )
    TrackBanner(track);

  // Pre tracking user intervention process.
  if( fpUserTrackingAction != 0 ) {
    fpUserTrackingAction->PreUserTrackingAction(track);
  }
//#ifdef G4_STORE_TRAJECTORY
  G4TrackingInformation* trackingInfo = GetIT(track)->GetTrackingInfo();
  G4Trajectory_Lock* trajectory_lock =
      trackingInfo->GetTrajectory_Lock();

  // Construct a trajectory if it is requested
  if(fStoreTrajectory&&(!trajectory_lock))
  {
    trajectory_lock = new G4Trajectory_Lock();
    trackingInfo->SetTrajectory_Lock(trajectory_lock);
    G4VTrajectory* trajectory = 0;
    // default trajectory concrete class object
    switch (fStoreTrajectory) {
      default:
      case 1: trajectory = new G4Trajectory(track); break;
      case 2: trajectory = new G4SmoothTrajectory(track); break;
      case 3: trajectory = new G4RichTrajectory(track); break;
    }
    trajectory_lock->fpTrajectory = trajectory;
  }
//#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
ITTrackingInteractivity::AppendStep(G4Track* track, G4Step* step)
{
  ////////////////////////////
  // If you want to use sensitive detector
  /*
    G4VPhysicalVolume* currentVolume =
    step->GetPreStepPoint()->GetPhysicalVolume();
    G4SteppingControl stepControlFlag =  step->GetControlFlag();

    if( currentVolume  != 0 && stepControlFlag != AvoidHitInvocation)
    {
        G4VSensitiveDetector* sensitive = step->GetPreStepPoint()->
                GetSensitiveDetector();
        if( sensitive != 0 )
        {
            sensitive->Hit(fpStep);
        }
    }
   */
  ////////////////////////////

  if(fpUserSteppingAction)
    fpUserSteppingAction->UserSteppingAction(step);

  ////////////////////////////
  // If you want to use regional stepping action
  /*
    G4UserSteppingAction* regionalAction
            = fpStep->GetPreStepPoint()->GetPhysicalVolume()->
                  GetLogicalVolume()->GetRegion()->
                  GetRegionalSteppingAction();
    if( regionalAction ) regionalAction->UserSteppingAction(fpStep);
   */
  ////////////////////////////

  if(fStoreTrajectory)
  {
    G4TrackingInformation* trackingInfo =
        GetIT(track)->GetTrackingInfo();
    G4Trajectory_Lock* trajectory_lock =
        trackingInfo->GetTrajectory_Lock();
    trajectory_lock->fpTrajectory->AppendStep(step);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ITTrackingInteractivity::EndTracking(G4Track* track)
{
#ifdef G4VERBOSE
  if(fVerboseLevel)
  {
    TrackBanner(track, "G4ITTrackingManager::EndTracking : ");
  }
#endif
  // Post tracking user intervention process.
  if( fpUserTrackingAction != 0 ) {
    fpUserTrackingAction->PostUserTrackingAction(track);
  }

//#ifdef G4_STORE_TRAJECTORY
  G4TrackingInformation* trackingInfo =
      GetIT(track)->GetTrackingInfo();
  G4Trajectory_Lock* trajectory_lock =
      trackingInfo->GetTrajectory_Lock();

  if(trajectory_lock)
  {
    G4VTrajectory*& trajectory = trajectory_lock->fpTrajectory;

    if(fStoreTrajectory && trajectory)
    {

#ifdef G4VERBOSE
      if(fVerboseLevel>10) trajectory->ShowTrajectory();
#endif
      G4TrackStatus istop = track->GetTrackStatus();

      if (trajectory && (istop != fStopButAlive) && (istop != fSuspend))
      {
        G4Event* currentEvent = G4EventManager::GetEventManager()
            ->GetNonconstCurrentEvent();

        if (currentEvent)
        {
          G4TrajectoryContainer* trajectoryContainer = currentEvent
              ->GetTrajectoryContainer();

          if (!trajectoryContainer)
          {
            trajectoryContainer = new G4TrajectoryContainer;
            currentEvent->SetTrajectoryContainer(trajectoryContainer);
          }
          trajectoryContainer->insert(trajectory);
        }
        else
        {
          fTrajectories.push_back(trajectory);
        }
      }
    }
    // Destruct the trajectory if it was created
    else if( (!fStoreTrajectory)&&trajectory ) {
      delete trajectory;
      trajectory = 0;
    }
    delete trajectory_lock;
    trackingInfo->SetTrajectory_Lock(0);
  }
//#endif
}

void ITTrackingInteractivity::Finalize()
{
  for (std::vector<G4VTrajectory*>::iterator it = fTrajectories.begin();
      it != fTrajectories.end(); it++)
  {
    G4VisManager::GetConcreteInstance()->Draw(**it);
  }
}
