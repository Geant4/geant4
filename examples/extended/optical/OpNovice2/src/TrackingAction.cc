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
/// \file optical/OpNovice2/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//

#include "TrackingAction.hh"
#include "TrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::TrackingAction()
:G4UserTrackingAction()
{;}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for track in tracking region
  TrackInformation* trackInfo = 
    (TrackInformation*)(aTrack->GetUserInformation());

  if (!trackInfo) {
    trackInfo = new TrackInformation(aTrack); 
    trackInfo->SetIsFirstTankX(true);
    aTrack->SetUserInformation(trackInfo);
  }

  trackInfo->SetIsFirstTankX(true);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if (secondaries)
  {
    TrackInformation* info = 
      (TrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if (nSeco > 0)
    {
      for(size_t i=0; i < nSeco; i++)
      { 
        TrackInformation* infoNew = new TrackInformation(info);
        (*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
}
