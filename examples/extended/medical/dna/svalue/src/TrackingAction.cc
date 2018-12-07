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
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"
#include "Run.hh"
#include "HistoManager.hh"
#include "MyFile.hh"

#ifdef MYFILE
 #include "MyPrimaryGeneratorActionFromFile.hh"
#else
 #include "PrimaryGeneratorAction.hh"
#endif

#include "G4RunManager.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifdef MYFILE

TrackingAction::TrackingAction(MyPrimaryGeneratorActionFromFile* prim)
:G4UserTrackingAction(),
 fPrimary(prim)
{}

#else

TrackingAction::TrackingAction(PrimaryGeneratorAction* prim)
:G4UserTrackingAction(),
 fPrimary(prim)
{}

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
 // histograms (useful for radionuclide spectrum, if set is vacuum)
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 if (aTrack->GetDefinition()==G4Electron::ElectronDefinition())
 {
   analysisManager->FillH1(0, aTrack->GetKineticEnergy());
   analysisManager->FillH1(1, aTrack->GetKineticEnergy());
 }
 if (aTrack->GetDefinition()==G4Gamma::GammaDefinition())
 {
   analysisManager->FillH1(2, aTrack->GetKineticEnergy());
   analysisManager->FillH1(3, aTrack->GetKineticEnergy());
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
 // histograms
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 
 G4int trackID = track->GetTrackID();
 
 Run* run
   = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
         
 //track length of primary particle or charged secondaries
 //
 G4double tracklen = track->GetTrackLength();
 if (trackID == 1) {
    run->AddTrackLength(tracklen);
    analysisManager->FillH1(7, tracklen);
 } else if (track->GetDefinition()->GetPDGCharge() != 0.)
    analysisManager->FillH1(10, tracklen);
           
 //extract projected range of primary particle
 //
 if (trackID == 1) {
   G4double pr = (track->GetPosition())*
                 (fPrimary->GetParticleGun()->GetParticleMomentumDirection());
   run->AddProjRange(pr);
   analysisManager->FillH1(9, pr);
 }
            
 //mean step size of primary particle
 //
 if (trackID == 1) {
   G4int nbOfSteps = track->GetCurrentStepNumber();
   G4double stepSize = tracklen/nbOfSteps;
   run->AddStepSize(nbOfSteps,stepSize);
 }
}
