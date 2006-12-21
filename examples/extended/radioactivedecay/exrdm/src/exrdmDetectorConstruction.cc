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
#include "exrdmDetectorConstruction.hh"
#include "exrdmDetectorMessenger.hh"
//#include "exrdmDetectorSD.hh"
#include "G4UImanager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
//#include "G4SDManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "exrdmMaterial.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
exrdmDetectorConstruction::exrdmDetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0), 
 solidDetector(0),logicDetector(0),physiDetector(0), 
 TargetMater(0), DetectorMater(0),
 fWorldLength(0.)
{
  detectorMessenger = new exrdmDetectorMessenger(this);
  DefineMaterials();
  fDetectorThickness = 2.* cm;
  fTargetRadius = 0.5 * cm;
  fDetectorLength = 5.0 * cm;      
  fTargetLength  = 1.0 * cm;         
//--------- Sizes of the principal geometrical components (solids)  ---------
  fWorldLength = std::max(fTargetLength,fDetectorLength);
  fWorldRadius = fTargetRadius + fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
exrdmDetectorConstruction::~exrdmDetectorConstruction()
{
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void exrdmDetectorConstruction::DefineMaterials()
{
//--------- Material definition ---------

  materialsManager = new exrdmMaterial();
  // Lead
  materialsManager->AddMaterial("Lead","Pb",11.3*g/cm3,"");
  //Germanium detector
  materialsManager->AddMaterial("Germanium","Ge",5.323*g/cm3,""); 
  //CsI
  materialsManager->AddMaterial("CsI","Cs-I",4.51*g/cm3,"");

  // G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    
  DefaultMater = materialsManager->GetMaterial("Air");
  TargetMater  = materialsManager->GetMaterial("CsI");
  DetectorMater = materialsManager->GetMaterial("Germanium");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* exrdmDetectorConstruction::Construct()
{
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  //--------- Sizes of the principal geometrical components (solids)  ---------

  fWorldLength = std::max(fTargetLength,fDetectorLength);
  fWorldRadius = fTargetRadius + fDetectorThickness;
   
  //------------------------------ 
  // World
  //------------------------------ 

 solidWorld= new G4Tubs("world",0.,fWorldRadius,fWorldLength/2.,0.,twopi);
 logicWorld= new G4LogicalVolume( solidWorld, DefaultMater, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
				 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = G4ThreeVector(0,0,0);
   
  solidTarget = new G4Tubs("target",0.,fTargetRadius,fTargetLength/2.,0.,twopi);
  logicTarget = new G4LogicalVolume(solidTarget,TargetMater,"Target",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  logicTarget,     // its logical volume				  
				  "Target",        // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //  G4cout << "Target is a cylinder with rdius of " << targetradius/cm << " cm of " 
  //       << TargetMater->GetName() << G4endl;

  //------------------------------ 
  // Detector
  //------------------------------
  
  G4ThreeVector positionDetector = G4ThreeVector(0,0,0);
  
  solidDetector = new G4Tubs("detector",fTargetRadius,fWorldRadius,fDetectorLength/2.,0.,twopi);
  logicDetector = new G4LogicalVolume(solidDetector ,DetectorMater, "Detector",0,0,0);  
  physiDetector = new G4PVPlacement(0,              // no rotation
				  positionDetector, // at (x,y,z)
				  logicDetector,    // its logical volume				  
				  "Detector",       // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  //  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // G4String detectortargetSDname = "exrdm/DetectorTargetSD";
  // exrdmDetectorSD* aDetectorSD = new exrdmDetectorSD( detectorTargetSDname );
  // SDman->AddNewDetector( aDetectorSD );
  //logicTarget->SetSensitiveDetector( aDetectorSD );
  // logicDetector->SetSensitiveDetector( aDetectorSD );
  //
  //-------------------------------------------------
  // regions
  //
  //  if(targetRegion) delete targetRegion;
  // if(detectorRegion) delete detectorRegion;
  targetRegion = new G4Region("Target");
  detectorRegion   = new G4Region("Detector");
  targetRegion->AddRootLogicalVolume(logicTarget);
  detectorRegion->AddRootLogicalVolume(logicDetector);

  //--------- Visualization attributes -------------------------------
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* TargetVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicTarget ->SetVisAttributes(TargetVisAtt);
  G4VisAttributes* DetectorVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,.0));
  logicDetector->SetVisAttributes(DetectorVisAtt);


  //------------ set the incident position ------

 // get the pointer to the User Interface manager 
    
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  //      UI->ApplyCommand("/run/verbose 1");
  //      UI->ApplyCommand("/event/verbose 2");
  //      UI->ApplyCommand("/tracking/verbose 1");  

  G4double zpos = -fWorldLength/2.;
  G4String command = "/gps/pos/centre ";
  std::ostringstream os;
  os << zpos ; 
  G4String xs = os.str();
  UI->ApplyCommand(command+"0. 0. "+xs+" mm");
  UI->ApplyCommand("/gps/pos/type Point");
  command = "/gps/position ";
  //  UI->ApplyCommand(command+"0. 0. "+xs+" mm");
  UI->ApplyCommand("/gps/particle proton");
  UI->ApplyCommand("/gps/direction 0 0 1");
  UI->ApplyCommand("/gps/energy 100 MeV");
  //       
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void exrdmDetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {TargetMater = pttoMaterial;
      if (logicTarget) logicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target has been changed to " << fTargetLength/cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmDetectorConstruction::setDetectorMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {DetectorMater = pttoMaterial;
      if (logicDetector) logicDetector->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The Deetctor has been changed to" << fDetectorLength/cm << " cm of "
             << materialName << G4endl;
     }             
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
