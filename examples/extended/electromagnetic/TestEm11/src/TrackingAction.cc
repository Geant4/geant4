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
// $Id: TrackingAction.cc,v 1.1 2005/06/03 15:20:32 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, RunAction* run,
                               HistoManager* histo)
:detector(det),runAction(run), histoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
 G4int trackID = track->GetTrackID();
  
 //track length of primary particle or charged secondaries
 //
 G4double tracklen = track->GetTrackLength();
 if (trackID == 1) {
    runAction->AddTrackLength(tracklen);
    histoManager->FillHisto(3, tracklen);
 } else if (track->GetDefinition()->GetPDGCharge() != 0.)
    histoManager->FillHisto(6, tracklen);
           
 //extract projected range of primary particle
 //
 if (trackID == 1) {
   G4double x = track->GetPosition().x() + 0.5*detector->GetAbsorSizeX();
   runAction->AddProjRange(x);
   histoManager->FillHisto(5, x);
 }
            
 //mean step size of primary particle
 //
 if (trackID == 1) {
   G4int nbOfSteps = track->GetCurrentStepNumber();
   G4double stepSize = tracklen/nbOfSteps;
   runAction->AddStepSize(nbOfSteps,stepSize);
 }
            
 //status of primary particle : absorbed, transmited, reflected ?
 //
 if (trackID == 1) {
  G4int status = 0;
  if (!track->GetNextVolume()) {
    if (track->GetMomentumDirection().x() > 0.) status = 1;
    else                                        status = 2;
  }
  runAction->AddTrackStatus(status);
 }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

