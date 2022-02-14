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
//  History :
//     21/10/2021 : DLa update to modify the material and the size
// -------------------------------------------------------------------
// -------------------------------------------------------------------

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(), fPhysiWorld(nullptr), fLogicWorld(nullptr),    fSolidWorld(nullptr),  fWorldMaterial(nullptr),  flogicTarget(nullptr),  ftargetSolid(nullptr),  fMaterial(nullptr),   fRegion(nullptr)
{
  fBoxSize = 1.0*um;
  DefineMaterials();
  SetMaterial("G4_Si");  
  fDetectorMessenger = new DetectorMessenger(this);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{
  delete fPhysiWorld;
  delete fLogicWorld;  
  delete fSolidWorld;  
  //delete fWorldMaterial;   // no delete because link to database
  delete flogicTarget;
  delete ftargetSolid;
  delete fMaterial; 
  delete fRegion;
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{ 

  // Silicon is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();
  fMaterial =  man->FindOrBuildMaterial("G4_Si");
  fWorldMaterial =  man->FindOrBuildMaterial("G4_Si");  // world material;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = 
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  
  if (pttoMaterial) {
    fMaterial = pttoMaterial;
    if ( flogicTarget ) { flogicTarget->SetMaterial(fMaterial); }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  // WORLD VOLUME
  
  
    G4double TargetSizeX =  1*um; 
  G4double TargetSizeY =  TargetSizeX; 
  G4double TargetSizeZ =  TargetSizeX; 
  fBoxSize = TargetSizeX;

  fWorldSizeX  = TargetSizeX*2; 
  fWorldSizeY  = TargetSizeY*2; 
  fWorldSizeZ  = TargetSizeZ*2; 

  fSolidWorld = new G4Box("World",				     //its name
			   fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2);  //its size
  

  fLogicWorld = new G4LogicalVolume(fSolidWorld,	//its solid
				   fWorldMaterial,		//its material
				   "World");		//its name
  
  fPhysiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 fLogicWorld,		//its logical volume
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number



  ftargetSolid = new G4Box("Target",				     //its name
				 TargetSizeX/2,TargetSizeY/2,TargetSizeZ/2);   //its size
  

  flogicTarget = new G4LogicalVolume(ftargetSolid,       //its solid
						     fMaterial,	//its material
						     "Target");		//its name
  
  new G4PVPlacement(0,			                               //no rotation
		    G4ThreeVector(),	                               //at (0,0,0)
		    "Target",		//its name
		    flogicTarget,	//its logical volume
		    fPhysiWorld,		//its mother  volume
		    false,		//no boolean operation
		    0);			//copy number

  // Visualization attributes
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  worldVisAtt->SetVisibility(true);
  fLogicWorld->SetVisAttributes(worldVisAtt);

  G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); 
  worldVisAtt1->SetVisibility(true);
  flogicTarget->SetVisAttributes(worldVisAtt1);

  // Create Target G4Region and add logical volume
  
  fRegion = new G4Region("Target");
  
  G4ProductionCuts* cuts = new G4ProductionCuts();
  
  G4double defCut = 1*nanometer;
  cuts->SetProductionCut(defCut,"gamma");
  cuts->SetProductionCut(defCut,"e-");
  cuts->SetProductionCut(defCut,"e+");
  cuts->SetProductionCut(defCut,"proton");
  
  fRegion->SetProductionCuts(cuts);
  fRegion->AddRootLogicalVolume(flogicTarget); 

  return fPhysiWorld;
}
 
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
  if(ftargetSolid) {
    ftargetSolid->SetXHalfLength(fBoxSize/2);
    ftargetSolid->SetYHalfLength(fBoxSize/2);
    ftargetSolid->SetZHalfLength(fBoxSize/2);
  }
  fWorldSizeX = value*2.0;
  fWorldSizeY = value*2.0;
  fWorldSizeZ = value*2.0;
  
  if(fSolidWorld) {
    fSolidWorld->SetXHalfLength(fWorldSizeX/2);
    fSolidWorld->SetYHalfLength(fWorldSizeY/2);
    fSolidWorld->SetZHalfLength(fWorldSizeZ/2);
  }
  
}

