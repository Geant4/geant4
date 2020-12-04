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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#include "SteppingAction.hh"
#include "Analysis.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* run,DetectorConstruction* det)
:fRun(run),fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* step)
  
{ 

G4AnalysisManager* man = G4AnalysisManager::Instance();


if  ( (step->GetTrack()->GetDynamicParticle()->GetDefinition() == 
       G4Proton::ProtonDefinition())

/*
// for doublet

    && (step->GetPostStepPoint()->GetPosition().z()/mm>-3230.2)
         && (step->GetPostStepPoint()->GetPosition().z()/mm<-3229.8) 
*/

// for triplet and whole line

         && (step->GetPostStepPoint()->GetPosition().z()/mm>249.99999)
         && (step->GetPostStepPoint()->GetPosition().z()/mm<250.00001) 
         && (step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->
                   GetLogicalVolume()->GetName()  == fDetector->GetLogicalVol()->GetName())
         && (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->
                  GetLogicalVolume()->GetName() == fDetector->GetLogicalWorld()->GetName())
    )
      
     {
         fXIn = step->GetPostStepPoint()->GetPosition().x();
         fYIn = step->GetPostStepPoint()->GetPosition().y();
         fZIn = step->GetPostStepPoint()->GetPosition().z();
         fE   = step->GetTrack()->GetKineticEnergy();

         G4ThreeVector angleIn;
         angleIn = step->GetTrack()->GetMomentumDirection();

         fThetaIn = std::asin(angleIn[0]/std::sqrt(angleIn[0]
                   *angleIn[0]+angleIn[1]*angleIn[1]+angleIn[2]*angleIn[2]));
         fPhiIn = std::asin(angleIn[1]/std::sqrt(angleIn[0]
                   *angleIn[0]+angleIn[1]*angleIn[1]+angleIn[2]*angleIn[2]));

         G4cout << "    =>IMAGE : X(microns)=" << fXIn/micrometer
                <<" Y(microns)="<< fYIn/micrometer << " THETA(mrad)="
                << (fThetaIn/mrad) << " PHI(mrad)=" << (fPhiIn/mrad) << G4endl;
         G4cout << G4endl;

         if (fDetector->GetCoef()==1) 
	 {
           fRun->AddRow();
           fRun->AddToXVector(fXIn/um);
           fRun->AddToYVector(fYIn/um);
           fRun->AddToThetaVector(fThetaIn/mrad);
           fRun->AddToPhiVector(fPhiIn/mrad);
         }

	 //Fill ntuple 3
	 man->FillNtupleDColumn(3,0,fXIn/um);
	 man->FillNtupleDColumn(3,1,fYIn/um);
	 man->FillNtupleDColumn(3,2,fThetaIn/mrad);
	 man->FillNtupleDColumn(3,3,fPhiIn/mrad);
	 man->AddNtupleRow(3);

     }

if (fDetector->GetProfile()==1) 
{

   if  (
            (step->GetTrack()->GetDynamicParticle()->GetDefinition()== G4Proton::ProtonDefinition())
         && (step->GetPreStepPoint()->GetTouchableHandle()
             ->GetVolume()->GetLogicalVolume()->GetName()  == fDetector->GetLogicalVol()->GetName())
         && (step->GetPostStepPoint()->GetTouchableHandle()
            ->GetVolume()->GetLogicalVolume()->GetName() == fDetector->GetLogicalVol()->GetName()) 
       )
   {
         fXIn = step->GetPostStepPoint()->GetPosition().x();
         fYIn = step->GetPostStepPoint()->GetPosition().y();
         fZIn = step->GetPostStepPoint()->GetPosition().z();

         //Fill ntuple 1
	 man->FillNtupleDColumn(1,0,fXIn/um);
	 man->FillNtupleDColumn(1,1,fYIn/um);
	 man->FillNtupleDColumn(1,2,fZIn/um);
	 man->AddNtupleRow(1);
   }
}
   
if (fDetector->GetGrid()==1) 
{

   if  (
            (step->GetTrack()->GetDynamicParticle()->GetDefinition()== G4Proton::ProtonDefinition())
         && (step->GetPreStepPoint()->GetTouchableHandle()
             ->GetVolume()->GetLogicalVolume()->GetName()  == fDetector->GetLogicalGrid()->GetName())
         && (step->GetPostStepPoint()->GetTouchableHandle()
             ->GetVolume()->GetLogicalVolume()->GetName() == fDetector->GetLogicalWorld()->GetName()) 
       )
   {
         fXIn = step->GetPostStepPoint()->GetPosition().x();
         fYIn = step->GetPostStepPoint()->GetPosition().y();
         fE   = step->GetTrack()->GetKineticEnergy();

	 //Fill ntuple 2
	 man->FillNtupleDColumn(2,0,fXIn/um);
	 man->FillNtupleDColumn(2,1,fYIn/um);
	 man->FillNtupleDColumn(2,2,fE/MeV);
	 man->AddNtupleRow(2);
   }
 }

// end
}     
