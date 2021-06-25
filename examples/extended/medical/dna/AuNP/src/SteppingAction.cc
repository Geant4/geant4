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
/// \file medical/dna/range/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 78723 2014-01-20 10:32:17Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4RunManager.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction(),
fpDetector(0),
fRNP(0),fRAbs(0)
{

  fpDetector =
      dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
                        ->GetUserDetectorConstruction());
  fRNP       = fpDetector->GetNPRadius()   /CLHEP::nm;
  fRAbs      = fpDetector->GetAbsRadius()  /CLHEP::nm;
  fTrackCut  = fpDetector->GetTrackingCut()/CLHEP::eV;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

 fpDetector =
     dynamic_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
                       ->GetUserDetectorConstruction());
 fRNP       = fpDetector->GetNPRadius()   /CLHEP::nm;
 fRAbs      = fpDetector->GetAbsRadius()  /CLHEP::nm;
 fTrackCut  = fpDetector->GetTrackingCut()/CLHEP::eV;

 G4ThreeVector pos     = aStep->GetPreStepPoint() ->GetPosition();
 G4ThreeVector postpos = aStep->GetPostStepPoint()->GetPosition();

 G4double R      = pos.mag  ()/CLHEP::nm;

 if(fRNP<R){
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   G4double energy = aStep->GetTotalEnergyDeposit()/CLHEP::joule;
   analysisManager->FillH1(1,R,energy);
 }

 R      = pos.mag()/CLHEP::nm;
 if(fRNP<R){
   if(std::abs(pos.z())<10*CLHEP::nm){
     G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
     G4double energy = aStep->GetTotalEnergyDeposit()/CLHEP::joule;
     G4double theta  = std::atan2(pos.y(),pos.x())/CLHEP::deg;
     if(0<=theta){
       analysisManager->FillH2(0,theta,R,energy);
     }else{
       analysisManager->FillH2(0,theta+360,R,energy);
     }
   }
 }
 
   if( aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() 
                                                   =="NanoParticle"){
      if( aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
        //*** WARNING: this line will kill all incident electrons 
        //*** at the end NP surface***
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        G4Track  *track  = aStep->GetTrack();
        G4double  trackE = track->GetKineticEnergy()/CLHEP::eV;

        if(aStep->GetTrack()->GetTrackID() == 1 ){
          if(pos.x()<0){
            analysisManager->FillH1(8,trackE);
          }else{
            analysisManager->FillH1(9,trackE);
          }
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }
        else{
          G4ThreeVector dir = aStep->GetPostStepPoint()->GetMomentumDirection();
          G4double dot      = dir.dot(postpos);
          if(dot<0.0){return;}

          if(trackE<fTrackCut){
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return;
          }

          if(track->GetDefinition()->GetPDGCharge() != 0){
            analysisManager->FillH1(4,trackE);
          }else {
            analysisManager->FillH1(5,trackE);
          }
        }
      }
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
