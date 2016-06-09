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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITTrackingManager.hh"
#include "G4Track.hh"
#include "G4ProcessManager.hh"
#include "G4IT.hh"
#include "G4TrackingInformation.hh"

/***
#include "G4UserTrackingAction.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4SmoothTrajectory.hh"
#include "G4RichTrajectory.hh"
***/

G4ITTrackingManager::G4ITTrackingManager()
{
/***
    G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
    fStoreTrajectory = trackingManager->GetStoreTrajectory();
    fVerboseLevel = trackingManager->GetVerboseLevel();
    fpUserTrackingAction = trackingManager->GetUserTrackingAction();
    fSetNewUserTrackingAction = false;
***/

    // if you uncomment the above part, please comment the part below
    fStoreTrajectory = 0;
    fVerboseLevel = 0;
    fSetNewUserTrackingAction = false;
}
//___________________________________________________
void G4ITTrackingManager::Initialize()
{
/***
    G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
    fStoreTrajectory = trackingManager->GetStoreTrajectory();
    fVerboseLevel = trackingManager->GetVerboseLevel();
    if(!fpUserTrackingAction)
        fpUserTrackingAction = trackingManager->GetUserTrackingAction();
***/
}

//___________________________________________________
/***
void G4ITTrackingManager::SetUserTrackingAction(G4UserTrackingAction* trackAct)
{
    G4TrackingManager* trackingManager = G4EventManager::GetEventManager()->GetTrackingManager();
    G4UserTrackingAction* std_trackAct = trackingManager->GetUserTrackingAction();

    if(trackAct != std_trackAct) fSetNewUserTrackingAction = true;
    else fSetNewUserTrackingAction = false;

    fpUserTrackingAction = trackAct;
}***/
//___________________________________________________
G4ITTrackingManager::~G4ITTrackingManager()
{
/***
    if(fSetNewUserTrackingAction)
        delete fpUserTrackingAction;
**/
}
//___________________________________________________
void G4ITTrackingManager::StartTracking(G4Track* track)
{
#ifdef G4VERBOSE
    if(fVerboseLevel)
    {
        TrackBanner(track, "G4ITTrackingManager::StartTracking : ");
    }
#endif

/***
    if(fVerboseLevel>0 && (G4VSteppingVerbose::GetSilent()!=1) ) TrackBanner(track);

    // Pre tracking user intervention process.
    if( fpUserTrackingAction != 0 ) {
        fpUserTrackingAction->PreUserTrackingAction(track);
    }
#ifdef G4_STORE_TRAJECTORY
    G4TrackingInformation* trackingInfo = GetIT(track)->GetTrackingInfo();
    G4Trajectory_Lock* trajectory_lock = trackingInfo->GetTrajectory_Lock();

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
#endif
***/

    // Inform beginning of tracking to physics processes
    track->GetDefinition()->GetProcessManager()->StartTracking(track);
}
//___________________________________________________
void G4ITTrackingManager::AppendTrajectory(G4Track* /***track***/, G4Step* /***step***/)
{
/***
#ifdef G4_STORE_TRAJECTORY
    if(fStoreTrajectory)
    {
        G4TrackingInformation* trackingInfo = GetIT(track)->GetTrackingInfo();
        G4Trajectory_Lock* trajectory_lock = trackingInfo->GetTrajectory_Lock();
        trajectory_lock->fpTrajectory->AppendStep(step);
    }
#endif
***/
}
//___________________________________________________
void G4ITTrackingManager::EndTracking(G4Track* track)
{
#ifdef G4VERBOSE
    if(fVerboseLevel)
    {
        TrackBanner(track, "G4ITTrackingManager::EndTracking : ");
    }
#endif
    // Post tracking user intervention process.
/***
    if( fpUserTrackingAction != 0 ) {
        fpUserTrackingAction->PostUserTrackingAction(track);
    }

#ifdef G4_STORE_TRAJECTORY
    G4TrackingInformation* trackingInfo = GetIT(track)->GetTrackingInfo();
    G4Trajectory_Lock* trajectory_lock = trackingInfo->GetTrajectory_Lock();

    if(trajectory_lock)
    {
        G4VTrajectory*& trajectory = trajectory_lock->fpTrajectory;

        if(fStoreTrajectory && trajectory)
        {

#ifdef G4VERBOSE
            if(fVerboseLevel>10) trajectory->ShowTrajectory();
#endif
            G4TrackStatus istop = track->GetTrackStatus();
            if(trajectory&&(istop!=fStopButAlive)&&(istop!=fSuspend))
            {
                G4Event* currentEvent = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
                G4TrajectoryContainer* trajectoryContainer = currentEvent->GetTrajectoryContainer();
                if(!trajectoryContainer)
                {
                    trajectoryContainer = new G4TrajectoryContainer;
                    currentEvent->SetTrajectoryContainer(trajectoryContainer);
                }
                trajectoryContainer->insert(trajectory);
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
#endif
***/
}


void G4ITTrackingManager::TrackBanner(G4Track* track, const G4String& message)
{
    G4cout << G4endl;
    G4cout << "*******************************************************"
           << "**************************************************"
           << G4endl;
    if(message != "")
        G4cout << message ;
    G4cout << " * G4Track Information: "
           << "   Particle : " << track->GetDefinition()->GetParticleName()
           << ","
           << "   Track ID : " << track->GetTrackID()
           << ","
           << "   Parent ID : " << track->GetParentID()
           << G4endl;
    G4cout << "*******************************************************"
           << "**************************************************"
           << G4endl;
    G4cout << G4endl;
}

