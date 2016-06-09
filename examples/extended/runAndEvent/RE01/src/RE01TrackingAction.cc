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
// $Id: RE01TrackingAction.cc,v 1.1 2004/11/26 07:37:43 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//



#include "RE01TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "RE01Trajectory.hh"
#include "RE01TrackInformation.hh"

void RE01TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for track in tracking region
  RE01TrackInformation* trackInfo = (RE01TrackInformation*)(aTrack->GetUserInformation());

  if(trackInfo->GetTrackingStatus() > 0)
  {
    fpTrackingManager->SetStoreTrajectory(true);
    fpTrackingManager->SetTrajectory(new RE01Trajectory(aTrack));
  }
  else
  { fpTrackingManager->SetStoreTrajectory(false); }
}

void RE01TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    RE01TrackInformation* info = (RE01TrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      { 
        RE01TrackInformation* infoNew = new RE01TrackInformation(info);
        (*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
}


