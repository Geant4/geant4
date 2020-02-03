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
/// \file runAndEvent/RE01/src/RE01TrackingAction.cc
/// \brief Implementation of the RE01TrackingAction class
//
//

#include "RE01TrackingAction.hh"
#include "RE01Trajectory.hh"
#include "RE01TrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01TrackingAction::RE01TrackingAction()
:G4UserTrackingAction()
{;}

void RE01TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for track in tracking region
  RE01TrackInformation* trackInfo = 
    (RE01TrackInformation*)(aTrack->GetUserInformation());

  if(trackInfo->GetTrackingStatus() > 0)
  {
    fpTrackingManager->SetStoreTrajectory(true);
    fpTrackingManager->SetTrajectory(new RE01Trajectory(aTrack));
  }
  else
  { fpTrackingManager->SetStoreTrajectory(false); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    RE01TrackInformation* info = 
      (RE01TrackInformation*)(aTrack->GetUserInformation());
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


