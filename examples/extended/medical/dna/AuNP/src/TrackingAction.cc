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
/// \file medical/dna/range/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id: TrackingAction.cc 78723 2014-01-20 10:32:17Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(PrimaryGeneratorAction* /*prim*/)
:G4UserTrackingAction(),
 //fPrimary(prim),
 fpDetector(0)
{
  fpDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
 fpDetector =
     dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
         ->GetUserDetectorConstruction());

 G4int trackID = track->GetTrackID();
 if (trackID == 1) fTrackLength      = 0;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

 G4ThreeVector pos     = track->GetPosition();
 G4double RNP = fpDetector->GetNPRadius()/CLHEP::nm;
 G4double R   = std::sqrt(pos.x()*pos.x()+pos.y()*pos.y()+pos.z()*pos.z())/CLHEP::nm;
 //G4double ene   = track->GetKineticEnergy();
 G4double trackE = track->GetKineticEnergy()/CLHEP::eV;
 if(RNP>R && track->GetTrackID() != 1 ){
   if(track->GetDefinition()->GetPDGCharge() != 0){
     analysisManager->FillH1(2,trackE);
   }else {
     analysisManager->FillH1(3,trackE);
   }
 }
 if(RNP<R && track->GetTrackID() != 1 ){
   if(track->GetDefinition()->GetPDGCharge() != 0){
     analysisManager->FillH1(6,R);
     analysisManager->FillH2(1,R,trackE);
   }else {
     analysisManager->FillH1(7,R);
     analysisManager->FillH2(2,R,trackE);
   }

 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* /*track*/)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
