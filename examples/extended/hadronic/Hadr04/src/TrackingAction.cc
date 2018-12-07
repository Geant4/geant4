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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction()
:G4UserTrackingAction(),
 fNbStep1(0),fNbStep2(0),fTrackLen1(0.),fTrackLen2(0.),fTime1(0.),fTime2(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track*)
{
  fNbStep1 = fNbStep2 = 0;
  fTrackLen1 = fTrackLen2 = 0.;
  fTime1 = fTime2 = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::UpdateTrackInfo(G4double ekin,G4double trackl,
                                     G4double time)
{
  const G4double thermal = 1*eV;
  if (ekin > thermal) {
    fNbStep1++; fTrackLen1 = trackl; fTime1 = time;    
  } else {
    fNbStep2++; fTrackLen2 = trackl - fTrackLen1; fTime2 = time - fTime1;  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
 // keep only primary neutron
 //
 G4int trackID = track->GetTrackID();
 if (trackID > 1) return;
 
 Run* run 
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
 run->SumTrackLength(fNbStep1,fNbStep2,fTrackLen1,fTrackLen2,fTime1,fTime2);
 
 // histograms
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->FillH1(1,fNbStep1);
 analysisManager->FillH1(2,fTrackLen1);
 analysisManager->FillH1(3,fTime1); 
 analysisManager->FillH1(4,fNbStep2);
 analysisManager->FillH1(5,fTrackLen2);
 analysisManager->FillH1(6,fTime2);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

