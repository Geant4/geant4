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
/// \file radioactivedecay/rdecay02/src/exrdmDetectorConstruction.cc
/// \brief Implementation of the exrdmDetectorConstruction class
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ios.hh"
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
exrdmDetectorConstruction::exrdmDetectorConstruction()
:fSolidWorld(0),  fLogicWorld(0),  fPhysiWorld(0),
 fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0), 
 fSolidDetector(0),fLogicDetector(0),fPhysiDetector(0),
 fMaterialsManager(0),
 fDefaultMater(0),fTargetMater(0),fDetectorMater(0),
 fTargetLength (1.*cm), fTargetRadius(0.5*cm),
 fDetectorLength(5.0 * cm), fDetectorThickness(2.0 * cm),
 fWorldLength (std::max(fTargetLength,fDetectorLength)),
 fWorldRadius (fTargetRadius + fDetectorThickness),
 fTargetRegion(0), fDetectorRegion(0)
{
  fDetectorMessenger = new exrdmDetectorMessenger(this);
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
exrdmDetectorConstruction::~exrdmDetectorConstruction()
{
  delete fDetectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void exrdmDetectorConstruction::DefineMaterials()
{
//--------- Material definition ---------

  fMaterialsManager = new exrdmMaterial();
  // Lead
  fMaterialsManager->AddMaterial("Lead","Pb",11.3*g/cm3,"");
  //Germanium detector
  fMaterialsManager->AddMaterial("Germanium","Ge",5.323*g/cm3,""); 
  //CsI
  fMaterialsManager->AddMaterial("CsI","Cs-I",4.51*g/cm3,"");

  // G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
    
  fDefaultMater = fMaterialsManager->GetMaterial("Air");
  fTargetMater  = fMaterialsManager->GetMaterial("CsI");
  fDetectorMater = fMaterialsManager->GetMaterial("Germanium");
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

 fSolidWorld= new G4Tubs("world",0.,fWorldRadius,fWorldLength/2.,0.,twopi);
 fLogicWorld= new G4LogicalVolume( fSolidWorld, fDefaultMater, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  fPhysiWorld = new G4PVPlacement(0,              // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 fLogicWorld,     // its logical volume
                                 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
                                 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = G4ThreeVector(0,0,0);
   
  fSolidTarget = new G4Tubs("target",0.,fTargetRadius,fTargetLength/2.,0.,twopi);
  fLogicTarget = new G4LogicalVolume(fSolidTarget,fTargetMater,"Target",0,0,0);
  fPhysiTarget = new G4PVPlacement(0,              // no rotation
                                  positionTarget,  // at (x,y,z)
                                  fLogicTarget,    // its logical volume                                  
                                  "Target",        // its name
                                  fLogicWorld,     // its mother  volume
                                  false,           // no boolean operations
                                  0);              // no particular field 

  //  G4cout << "Target is a cylinder with rdius of " << targetradius/cm << " cm of " 
  //       << fTargetMater->GetName() << G4endl;

  //------------------------------ 
  // Detector
  //------------------------------
  
  G4ThreeVector positionDetector = G4ThreeVector(0,0,0);
  
  fSolidDetector = new G4Tubs("detector",fTargetRadius,fWorldRadius,
                                                          fDetectorLength/2.,0.,twopi);
  fLogicDetector = new G4LogicalVolume(fSolidDetector ,fDetectorMater,
                                                                      "Detector",0,0,0);
  fPhysiDetector = new G4PVPlacement(0,             // no rotation
                                  positionDetector, // at (x,y,z)
                                  fLogicDetector,   // its logical volume                                  
                                  "Detector",       // its name
                                  fLogicWorld,      // its mother  volume
                                  false,            // no boolean operations
                                  0);               // no particular field 

  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  //  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // G4String detectortargetSDname = "exrdm/DetectorTargetSD";
  // exrdmDetectorSD* aDetectorSD = new exrdmDetectorSD( detectorTargetSDname );
  // SDman->AddNewDetector( aDetectorSD );
  //fLogicTarget->SetSensitiveDetector( aDetectorSD );
  // fLogicDetector->SetSensitiveDetector( aDetectorSD );
  //
  //-------------------------------------------------
  // regions
  //
  //  if(fTargetRegion) delete fTargetRegion;
  // if(fDetectorRegion) delete fDetectorRegion;
  fTargetRegion = new G4Region("Target");
  fDetectorRegion   = new G4Region("Detector");
  fTargetRegion->AddRootLogicalVolume(fLogicTarget);
  fDetectorRegion->AddRootLogicalVolume(fLogicDetector);

  //--------- Visualization attributes -------------------------------
  fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* TargetVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fLogicTarget ->SetVisAttributes(TargetVisAtt);
  G4VisAttributes* DetectorVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,.0));
  fLogicDetector->SetVisAttributes(DetectorVisAtt);


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
  
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void exrdmDetectorConstruction::SetTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {fTargetMater = pttoMaterial;
      if (fLogicTarget) fLogicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target has been changed to " << fTargetLength/cm
                     << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmDetectorConstruction::SetDetectorMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {fDetectorMater = pttoMaterial;
      if (fLogicDetector) fLogicDetector->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The Deetctor has been changed to" << fDetectorLength/cm
                     << " cm of "
             << materialName << G4endl;
     }             
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
