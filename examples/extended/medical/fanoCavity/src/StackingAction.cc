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
/// \file medical/fanoCavity/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4RunManager.hh"
#include "StackingMessenger.hh"
#include "Run.hh"

#include "G4Track.hh"
#include "G4EmCalculator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(DetectorConstruction* det)
:fDetector(det),fStackMessenger(0)
{
  fMatWall = 0;
  fZcav = 0.; 
  fEmCal = 0;
  first = true;
  fKillTrack  = true;
  
  //create a messenger for this class  
  fStackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete fEmCal;
  delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
 //get fDetector informations
 if (first) {
   fMatWall = fDetector->GetWallMaterial();
   fZcav    = 0.5*(fDetector->GetCavityThickness());
   fEmCal   = new G4EmCalculator();
   first   = false;
 }

 G4ClassificationOfNewTrack status = fUrgent;  

 //keep primary particle or neutral
 //
 G4ParticleDefinition* particle = track->GetDefinition();
 G4bool neutral = (particle->GetPDGCharge() == 0.);
 if ((track->GetParentID() == 0) || neutral) return status;

 Run* run = static_cast<Run*>(
              G4RunManager::GetRunManager()->GetNonConstCurrentRun());

 //energy spectrum of charged secondaries
 //
  G4double energy = track->GetKineticEnergy();  
  run->SumEsecond(energy);
  
 // kill e- which cannot reach Cavity
 //
 G4double position = (track->GetPosition()).z();
 G4double safe = std::abs(position) - fZcav;
 G4double range = fEmCal->GetRangeFromRestricteDEDX(energy,particle,fMatWall);
 if (fKillTrack) {
   if (range < 0.8*safe) status = fKill;
 }

 //histograms
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
 analysisManager->FillH1(1,position);
 analysisManager->FillH1(2,energy);
 G4ThreeVector direction = track->GetMomentumDirection();
 analysisManager->FillH1(3,std::acos(direction.z()));      
    
 return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
