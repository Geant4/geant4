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
// $Id: TrackingAction.cc,v 1.2 2006-06-29 16:41:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

