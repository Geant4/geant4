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
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "Analysis.hh"

#include "G4Alpha.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(RunAction* run,DetectorConstruction* det)
:fRun(run),fDetector(det)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 
  // Analysis manager
  
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Read phantom - Singleton
   
  fMyCellParameterisation = CellParameterisation::Instance(); 

  //

  // Material : 1 is cytoplasm, 2 is nucleus

  G4int matVoxelPRE = -1;
  G4int matVoxelPOST = -1;
  G4int tmp=-1;
  
  tmp = aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(); 
    
  if (tmp>0)  
  {
    matVoxelPRE =  fMyCellParameterisation->GetTissueType(tmp);
  }
  
  tmp = aStep->GetPostStepPoint()->GetTouchableHandle()->GetReplicaNumber();
  
  if (tmp>0) 
  {
    matVoxelPOST =  fMyCellParameterisation->GetTissueType(tmp);
  }
  
// COUNT GAS DETECTOR HITS

if (       ((aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalCollDetYoke())
        &&  (aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalIsobutane())
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition() ))
	
        || 
	   ((aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalCollDetGap4())
        &&  (aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalIsobutane())
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition() ))

	||

    	   ((aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalCollDetGap4())
        &&  (aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalIsobutane())
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition() ))

    )
	{
	  fRun->AddNbOfHitsGas();	
	}
	
// STOPPING POWER AND BEAM SPOT SIZE AT CELL ENTRANCE

if (       ((aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalPolyprop())
        &&  (aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalKgm())
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition() ))
	
        || 
	   ((aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalPolyprop())
        &&  (matVoxelPOST == 1)
        &&  (aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition() ))
    )
	{
	
	 if( (aStep->GetPreStepPoint()->GetKineticEnergy() - aStep->GetPostStepPoint()->GetKineticEnergy() ) >0) 
	 {
	   //Fill ntupleid=1 
	   man->FillNtupleDColumn(1,0,aStep->GetPreStepPoint()->GetKineticEnergy()/keV);
	   man->FillNtupleDColumn(1,1,
				  (aStep->GetPreStepPoint()->GetKineticEnergy() -
				   aStep->GetPostStepPoint()->GetKineticEnergy())/
				  keV/(aStep->GetStepLength()/micrometer));
	   man->AddNtupleRow(1);
	 }

         // Average dE over step suggested by Michel Maire

	 G4StepPoint* p1 = aStep->GetPreStepPoint();
         G4ThreeVector coord1 = p1->GetPosition();
         const G4AffineTransform transformation1 = p1->GetTouchable()->GetHistory()->GetTopTransform();
         G4ThreeVector localPosition1 = transformation1.TransformPoint(coord1);

	 G4StepPoint* p2 = aStep->GetPostStepPoint();
         G4ThreeVector coord2 = p2->GetPosition();
         const G4AffineTransform transformation2 = p2->GetTouchable()->GetHistory()->GetTopTransform();
         G4ThreeVector localPosition2 = transformation2.TransformPoint(coord2);

         G4ThreeVector localPosition = localPosition1 + G4UniformRand()*(localPosition2-localPosition1);
	 
	 // end
	 
	 //Fill ntupleid=2
	 man->FillNtupleDColumn(2,0,localPosition.x()/micrometer);
	 man->FillNtupleDColumn(2,1,localPosition.y()/micrometer);
	 man->AddNtupleRow(2);
	}

// ALPHA RANGE

if (

	(aStep->GetTrack()->GetDynamicParticle()->GetDefinition() == G4Alpha::AlphaDefinition())
	
	&&
	
	(aStep->GetTrack()->GetKineticEnergy()<1e-6)
	
	&&
			
          ( (matVoxelPOST==1)
	||  (aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume() == fDetector->GetLogicalKgm())	
	||  (matVoxelPOST==2) )	
		
   )
	
	{
	   //Fill ntupleid=3
	  man->FillNtupleDColumn(3,0,
				 aStep->GetPostStepPoint()->GetPosition().x()/micrometer);
	  man->FillNtupleDColumn(3,1,
				 aStep->GetPostStepPoint()->GetPosition().y()/micrometer);
	  man->FillNtupleDColumn(3,2,
				 aStep->GetPostStepPoint()->GetPosition().z()/micrometer);
	  man->AddNtupleRow(3);
 	}

// TOTAL DOSE DEPOSIT AND DOSE DEPOSIT WITHIN A PHANTOM VOXEL
// FOR ALL PARTICLES

if (matVoxelPRE  == 2)

	{ 
   	 G4double dose = (aStep->GetTotalEnergyDeposit()/joule)/(fRun->GetMassNucleus()/kg);
   	 fRun->AddDoseN(dose);

	 G4ThreeVector v;
    	 fRun->AddDoseBox(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(),
	  aStep->GetTotalEnergyDeposit()/eV);
	}


if (matVoxelPRE  == 1)

	{ 
   	 G4double dose = (aStep->GetTotalEnergyDeposit()/joule)/(fRun->GetMassCytoplasm()/kg);
   	 fRun->AddDoseC(dose);

	 G4ThreeVector v;
    	 fRun->AddDoseBox(aStep->GetPreStepPoint()->GetTouchableHandle()->GetReplicaNumber(),
	  aStep->GetTotalEnergyDeposit()/eV);
 	}

// PROTECTION AGAINST MSC LOOPS FOR e-

if ( aStep->GetTotalEnergyDeposit()/MeV<1e-25
     && aStep->GetTrack()->GetDefinition()==G4Electron::ElectronDefinition()) 
{
  aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  /*
  G4cout << "*** Warning *** : msc loop for " 
    << aStep->GetTrack()->GetDefinition()->GetParticleName() 
    << " in " << 
    aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
  */
}

}
