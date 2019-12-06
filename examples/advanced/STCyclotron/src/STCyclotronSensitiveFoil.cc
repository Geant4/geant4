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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronSensitiveFoil.cc
//
#include "STCyclotronRun.hh"
#include "STCyclotronSensitiveFoil.hh"
#include "STCyclotronAnalysis.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

STCyclotronSensitiveFoil::STCyclotronSensitiveFoil(const G4String& name,
					           STCyclotronDetectorConstruction* det)
  : G4VSensitiveDetector(name), 
    fDet(det)
{ 
  fTempTrack = 0;
  fTempTrack1 = 0;
  fTempEnergy = 0.;
  fTempVector = G4ThreeVector(0.,0.,0.);  
  fRun =0; 
}

STCyclotronSensitiveFoil::~STCyclotronSensitiveFoil()
{ 
delete fRun;
}

G4bool STCyclotronSensitiveFoil::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  fRun  = static_cast<STCyclotronRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4Track* fTrack = aStep->GetTrack();

  auto analysisManager = G4AnalysisManager::Instance();
  
  //Step/track information
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ThreeVector momentumDirection = aStep->GetPreStepPoint()->GetMomentumDirection();
  G4ThreeVector vectorPosition = aStep->GetPreStepPoint()->GetPosition();
  G4String name = fTrack->GetDefinition()->GetParticleName();
  
  //Collect general information concerning all of the particles
  fRun->EnergyDepositionFoil(edep);

  //Collect information about protons
  if(name == "proton" || name == "deuteron"){

    if(fTrack->GetTrackID()!=fTempTrack && (momentumDirection.getZ()>0.) &&
       vectorPosition.getX()< 7.5 &&  
       vectorPosition.getX()>-7.5 &&
       vectorPosition.getY()< 7.5 &&  
       vectorPosition.getY()>-7.5){
    
      analysisManager->FillH2(1,vectorPosition.getX(),vectorPosition.getY());
      analysisManager->FillH1(1,energy);
    
      fTempTrack = fTrack->GetTrackID();
    
    }
    
    if(fTempTrack1 == 0){
      fTempTrack1 = fTrack->GetTrackID();
    }

    if(fTrack->GetTrackID()!=fTempTrack1 && (momentumDirection.getZ()>0.) &&
       vectorPosition.getX()< 7.5 &&  
       vectorPosition.getX()>-7.5 &&
       vectorPosition.getY()< 7.5 &&  
       vectorPosition.getY()>-7.5 ){

      analysisManager->FillH2(5,fTempVector.getX(),fTempVector.getY());
      analysisManager->FillH1(3,fTempEnergy);

      fTempTrack1 = fTrack->GetTrackID();

    }
 
    fTempVector =  aStep->GetPostStepPoint()->GetPosition();//vectorPosition;
    fTempEnergy = aStep->GetPostStepPoint()->GetKineticEnergy();//energy;
  }
  
  fRun->SetFoilVolume(fDet->GetFoilVolume());
  fRun->SetFoilThickness(fDet->GetFoilThickness());
  return true;
}
