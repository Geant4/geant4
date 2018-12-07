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
/// \file persistency/P01/src/ExP01DetectorConstruction.cc
/// \brief Implementation of the ExP01DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "ExP01DetectorConstruction.hh"
#include "ExP01DetectorMessenger.hh"
#include "ExP01ChamberParameterisation.hh"
#include "ExP01MagneticField.hh"
#include "ExP01TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExP01DetectorConstruction::ExP01DetectorConstruction()
:G4VUserDetectorConstruction(),
 fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
 fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0), 
 fSolidTracker(0), fLogicTracker(0), fPhysiTracker(0), 
 fSolidChamber(0), fLogicChamber(0), fPhysiChamber(0), 
 fTargetMater(0), fChamberMater(0), fPMagField(0), fDetectorMessenger(0),
 fWorldLength(0.), fTargetLength(0.), fTrackerLength(0.),
 fNbOfChambers(0), fChamberWidth(0.), fChamberSpacing(0.)
{
  fPMagField = new ExP01MagneticField();
  fDetectorMessenger = new ExP01DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExP01DetectorConstruction::~ExP01DetectorConstruction()
{
  delete fPMagField;
  delete fDetectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* ExP01DetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel;

  //Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  //Lead
  G4Material* Pb = 
  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
    
  //Xenon gas
  G4Material* Xenon = 
  new G4Material("XenonGas", z=54., a=131.29*g/mole, density= 5.458*mg/cm3,
   kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  fNbOfChambers = 5;
  fChamberWidth = 20*cm;
  fChamberSpacing = 80*cm;
  
  fTrackerLength = (fNbOfChambers+1)*fChamberSpacing; // Full length of Tracker
  fTargetLength  = 5.0 * cm;                        // Full length of Target
  
  fTargetMater  = Pb;
  fChamberMater = Xenon;
  
  fWorldLength= 1.2 *(fTargetLength+fTrackerLength);
   
  G4double targetSize  = 0.5*fTargetLength;    // Half length of the Target  
  G4double trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker
      
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 

  G4double HalfWorldLength = 0.5*fWorldLength;
  
 fSolidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
 fLogicWorld= new G4LogicalVolume( fSolidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  fPhysiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 fLogicWorld,      // its logical volume
                                 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // copy number
                                 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = G4ThreeVector(0,0,-(targetSize+trackerSize));
   
  fSolidTarget = new G4Box("target",targetSize,targetSize,targetSize);
  fLogicTarget = new G4LogicalVolume(fSolidTarget,fTargetMater,"Target",0,0,0);
  fPhysiTarget = new G4PVPlacement(0,               // no rotation
                                  positionTarget,  // at (x,y,z)
                                  fLogicTarget,     // its logical volume
                                  "Target",        // its name
                                  fLogicWorld,      // its mother  volume
                                  false,           // no boolean operations
                                  0);              // copy number 

  G4cout << "Target is " << fTargetLength/cm << " cm of " 
         << fTargetMater->GetName() << G4endl;

  //------------------------------ 
  // Tracker
  //------------------------------
  
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);
  
  fSolidTracker = new G4Box("tracker",trackerSize,trackerSize,trackerSize);
  fLogicTracker = new G4LogicalVolume(fSolidTracker , Air, "Tracker",0,0,0);  
  fPhysiTracker = new G4PVPlacement(0,              // no rotation
                                  positionTracker, // at (x,y,z)
                                  fLogicTracker,    // its logical volume
                                  "Tracker",       // its name
                                  fLogicWorld,      // its mother  volume
                                  false,           // no boolean operations
                                  0);              // copy number 

  //------------------------------ 
  // Tracker segments
  //------------------------------
  // 
  // An example of Parameterised volumes
  // dummy values for G4Box -- modified by parameterised volume

  fSolidChamber = new G4Box("chamber", 100*cm, 100*cm, 10*cm); 
  fLogicChamber = new G4LogicalVolume(fSolidChamber, fChamberMater,"Chamber",0,0,0);
  
  G4double firstPosition = -trackerSize + 0.5*fChamberWidth;
  G4double firstLength = fTrackerLength/10;
  G4double lastLength  = fTrackerLength;
   
  G4VPVParameterisation* chamberParam = new ExP01ChamberParameterisation(  
                           fNbOfChambers,          // NoChambers 
                           firstPosition,         // Z of center of first 
                           fChamberSpacing,        // Z spacing of centers
                           fChamberWidth,          // Width Chamber 
                           firstLength,           // lengthInitial 
                           lastLength);           // lengthFinal
                           
  // dummy value : kZAxis -- modified by parameterised volume
  //
  fPhysiChamber = new G4PVParameterised(
                            "Chamber",       // their name
                            fLogicChamber,    // their logical volume
                            fLogicTracker,    // Mother logical volume
                            kZAxis,          // Are placed along this axis 
                            fNbOfChambers,    // Number of chambers
                            chamberParam);   // The parametrisation

  G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
         << "The chambers are " << fChamberWidth/mm << " mm of " 
         << fChamberMater->GetName() << "\n The distance between chamber is "
         << fChamberSpacing/cm << " cm" << G4endl;
         
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "ExP01/TrackerChamberSD";
  ExP01TrackerSD* aTrackerSD = new ExP01TrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  fLogicChamber->SetSensitiveDetector( aTrackerSD );

//--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  fLogicWorld  ->SetVisAttributes(BoxVisAtt);  
  fLogicTarget ->SetVisAttributes(BoxVisAtt);
  fLogicTracker->SetVisAttributes(BoxVisAtt);
  
  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  fLogicChamber->SetVisAttributes(ChamberVisAtt);
  
//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in N02PhysicsList how to setup the processes
  // G4StepLimiter or G4UserSpecialCuts).
    
  // Sets a max Step length in the tracker region, with G4StepLimiter
  //
  G4double maxStep = 0.5*fChamberWidth; 
  fLogicTracker->SetUserLimits(new G4UserLimits(maxStep));
  
  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExP01DetectorConstruction::SetTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {fTargetMater = pttoMaterial;
      fLogicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExP01DetectorConstruction::SetChamberMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {fChamberMater = pttoMaterial;
      fLogicChamber->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The chambers are " << fChamberWidth/cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExP01DetectorConstruction::SetMagField(G4double fieldValue)
{
  fPMagField->SetFieldValue(fieldValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
