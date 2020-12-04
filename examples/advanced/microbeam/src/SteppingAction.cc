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

SteppingAction::SteppingAction(RunAction* run,const DetectorConstruction* det)
:fRun(run),fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SteppingAction::UserSteppingAction(const G4Step* aStep)
  
{ 
  // Analysis manager
  
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  // Read phantom - Singleton
   
  fMyCellParameterisation = CellParameterisation::Instance(); 

  // Material : 1 is cytoplasm, 2 is nucleus

  G4int matVoxelPRE = -1;
  G4int matVoxelPOST = -1;

  const G4StepPoint* preStep = aStep->GetPreStepPoint();
  const G4StepPoint* postStep = aStep->GetPostStepPoint();
  const G4Track* track = aStep->GetTrack();

  const G4LogicalVolume* preVolume =
    preStep->GetPhysicalVolume()->GetLogicalVolume();

  const G4LogicalVolume* postVolume = nullptr;
  if(postStep->GetPhysicalVolume()) 
  {
    postVolume = postStep->GetPhysicalVolume()->GetLogicalVolume();
  }
  const G4ParticleDefinition* particle = 
    track->GetDynamicParticle()->GetDefinition();

  G4int preReplicaNumber = preStep->GetTouchableHandle()->GetReplicaNumber(); 
  G4double edep = aStep->GetTotalEnergyDeposit();
    
  if (preReplicaNumber>0)  
  {
    matVoxelPRE =  fMyCellParameterisation->GetTissueType(preReplicaNumber);
  }

  if(postVolume) 
  {
    G4int postReplicaNumber = postStep->GetTouchableHandle()->GetReplicaNumber();
    if (postReplicaNumber>0) 
    {
      matVoxelPOST = fMyCellParameterisation->GetTissueType(postReplicaNumber);
    }
  }
  
  // COUNT GAS DETECTOR HITS

  if (particle == G4Alpha::AlphaDefinition()) 
  { 
    if(postVolume == fDetector->GetLogicalIsobutane() &&
       ((preVolume == fDetector->GetLogicalCollDetYoke())
	|| 
	(preVolume == fDetector->GetLogicalCollDetGap4())
	||
	(preVolume == fDetector->GetLogicalCollDetGap4())))
    {
      fRun->AddNbOfHitsGas();	
    }
  
    // STOPPING POWER AND BEAM SPOT SIZE AT CELL ENTRANCE
    if(preVolume == fDetector->GetLogicalPolyprop() &&
       ( (postVolume == fDetector->GetLogicalKgm()) || 
	 (matVoxelPOST == 1)) )
    {	
      G4double deltaE = preStep->GetKineticEnergy() 
	- postStep->GetKineticEnergy();
      if(deltaE > 0.0) 
      {
	//Fill ntupleid=1 
	man->FillNtupleDColumn(1,0,preStep->GetKineticEnergy()/keV);
	man->FillNtupleDColumn(1,1,deltaE*micrometer/(keV*aStep->GetStepLength()));
	man->AddNtupleRow(1);
      }

      // Average dE over step suggested by Michel Maire
      G4ThreeVector coord1 = preStep->GetPosition();
      const G4AffineTransform transformation1 = 
	preStep->GetTouchable()->GetHistory()->GetTopTransform();
      G4ThreeVector localPosition1 = transformation1.TransformPoint(coord1);

      G4ThreeVector coord2 = postStep->GetPosition();
      const G4AffineTransform transformation2 = 
	postStep->GetTouchable()->GetHistory()->GetTopTransform();
      G4ThreeVector localPosition2 = transformation2.TransformPoint(coord2);

      G4ThreeVector localPosition = 
	localPosition1 + G4UniformRand()*(localPosition2-localPosition1);
	 
      //Fill ntupleid=2
      man->FillNtupleDColumn(2,0,localPosition.x()/micrometer);
      man->FillNtupleDColumn(2,1,localPosition.y()/micrometer);
      man->AddNtupleRow(2);
    }

    // ALPHA RANGE
    if (postStep->GetKineticEnergy() < eV && 
	( (matVoxelPOST==1) ||
	  (postVolume == fDetector->GetLogicalKgm()) ||	
	  (matVoxelPOST==2) ) )
    {
      //Fill ntupleid=3
      man->FillNtupleDColumn(3,0,postStep->GetPosition().x()/micrometer);
      man->FillNtupleDColumn(3,1,postStep->GetPosition().y()/micrometer);
      man->FillNtupleDColumn(3,2,postStep->GetPosition().z()/micrometer);
      man->AddNtupleRow(3);
    }
    
    // TOTAL DOSE DEPOSIT AND DOSE DEPOSIT WITHIN A PHANTOM VOXEL
    // FOR ALL PARTICLES
  }

  if (matVoxelPRE  == 2)
  { 
    G4double dose = (edep/joule)/(fRun->GetMassNucleus()/kg);
    fRun->AddDoseN(dose);
    fRun->AddDoseBox(preReplicaNumber, edep/eV);
  }
  else if (matVoxelPRE  == 1)
  { 
    G4double dose = (edep/joule)/(fRun->GetMassCytoplasm()/kg);
    fRun->AddDoseC(dose);
    fRun->AddDoseBox(preReplicaNumber, edep/eV);
  }

  // PROTECTION AGAINST POSSIBLE MSC LOOPS FOR e-

  //  if ( edep/MeV<1e-25 && particle == G4Electron::Electron()) 
  //  {
      //aStep->GetTrack()->SetTrackStatus(fStopAndKill);
      /*
	G4cout << "*** Warning *** : msc loop for " 
	<< track->GetDefinition()->GetParticleName() 
	<< " in " << 
	postPoint->GetTouchableHandle()->GetVolume()->GetName() << G4endl;
      */
  //  }
}
