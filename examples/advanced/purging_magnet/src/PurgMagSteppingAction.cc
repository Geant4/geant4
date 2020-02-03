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
// Code developed by:
//  S.Larsson
//
//    ***********************************
//    *                                 *
//    *    PurgMagSteppingAction.cc     *
//    *                                 *
//    ***********************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "PurgMagSteppingAction.hh"
#include "PurgMagRunAction.hh"
#include "PurgMagDetectorConstruction.hh"

#include "G4SteppingManager.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

#include "PurgMagAnalysisManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagSteppingAction::PurgMagSteppingAction(const 
					     PurgMagDetectorConstruction* det)
:Detector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PurgMagSteppingAction::~PurgMagSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PurgMagSteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  
  //Collection at SSD in N-tuples. Electrons and photons separated
  //Prestep point in World, next volume MeasureVolume, process transportation
  if ((aStep->GetPreStepPoint()->GetPhysicalVolume() == Detector->GetWorld())&&
      (aStep->GetTrack()->GetNextVolume() == Detector->GetMeasureVolume())&&
      //(aStep->GetTrack()->GetMomentumDirection().z()>0.)&& // only particles with positive momentum
      (aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "Transportation"))
    {
      G4double gx, gy, gz, ge, gpx, gpy, gpz, ex, ey, ez, ee; 
      G4double epx, epy, epz, px, py, pz, pe, ppx, ppy, ppz;   
		 
      // Electrons 
      if(aStep->GetTrack()->GetDynamicParticle()->GetDefinition()
	 == G4Electron::Definition())
	{//Position 
         ex = (aStep->GetTrack()->GetPosition().x())/cm;
         ey = (aStep->GetTrack()->GetPosition().y())/cm; 
         ez = (aStep->GetTrack()->GetPosition().z())/cm;
	 // Energy
	 ee = (aStep->GetTrack()->GetKineticEnergy())/MeV;
	 // Momentum
	 epx = aStep->GetTrack()->GetMomentum().x();
	 epy = aStep->GetTrack()->GetMomentum().y();
	 epz = aStep->GetTrack()->GetMomentum().z();
	  
	 // Fill N-tuple electrons  ( id = 1)
	 man->FillNtupleDColumn(1,0, ex);
	 man->FillNtupleDColumn(1,1, ey);
	 man->FillNtupleDColumn(1,2, ez);
	 man->FillNtupleDColumn(1,3, ee);
	 man->FillNtupleDColumn(1,4, epx);
	 man->FillNtupleDColumn(1,5, epy);
	 man->FillNtupleDColumn(1,6, epz);
	 man->AddNtupleRow(1);  
      }

      // Photons
      if (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == 
	  G4Gamma::Definition())
	{
	  
	  // Position
          gx = (aStep->GetTrack()->GetPosition().x())/cm;
	  gy = (aStep->GetTrack()->GetPosition().y())/cm;
	  gz = (aStep->GetTrack()->GetPosition().z())/cm;
	  
	  // Energy
	  ge = (aStep->GetTrack()->GetKineticEnergy())/MeV;
	  
	  // Momentum
	  gpx = aStep->GetTrack()->GetMomentum().x();
	  gpy = aStep->GetTrack()->GetMomentum().y();
	  gpz = aStep->GetTrack()->GetMomentum().z();

	  // Fill N-tuple photons (id=2)
	  man->FillNtupleDColumn(2,0, gx);
	  man->FillNtupleDColumn(2,1, gy);
	  man->FillNtupleDColumn(2,2, gz);
	  man->FillNtupleDColumn(2,3, ge);
	  man->FillNtupleDColumn(2,4, gpx);
	  man->FillNtupleDColumn(2,5, gpy);
	  man->FillNtupleDColumn(2,6, gpz);
	  man->AddNtupleRow(2); 
	}


      // Positrons
      if (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == 
	  G4Positron::Definition())
	{

	  // Position
          px = (aStep->GetTrack()->GetPosition().x())/cm;
	  py = (aStep->GetTrack()->GetPosition().y())/cm;
	  pz = (aStep->GetTrack()->GetPosition().z())/cm;
	  
	  // Energy
	  pe = (aStep->GetTrack()->GetKineticEnergy())/MeV;
	  
	  // Momentum
	  ppx = aStep->GetTrack()->GetMomentum().x();
	  ppy = aStep->GetTrack()->GetMomentum().y();
	  ppz = aStep->GetTrack()->GetMomentum().z();
	  
	  // Fill Ntuple positrons ( id = 3)
	  man->FillNtupleDColumn(3,0, px);
	  man->FillNtupleDColumn(3,1, py);
	  man->FillNtupleDColumn(3,2, pz);
	  man->FillNtupleDColumn(3,3, pe);
	  man->FillNtupleDColumn(3,4, ppx);
	  man->FillNtupleDColumn(3,5, ppy);
	  man->FillNtupleDColumn(3,6, ppz);
	  man->AddNtupleRow(3);  
	}
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
