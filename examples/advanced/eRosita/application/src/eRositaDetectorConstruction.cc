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
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "eRositaDetectorConstruction.hh"
#include "eRositaTrackerSD.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
eRositaDetectorConstruction::eRositaDetectorConstruction()

{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
eRositaDetectorConstruction::~eRositaDetectorConstruction()
{
             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* eRositaDetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, z;
  G4double density;
//   G4int nel;

//   //Air
//   G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
//   G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
//   G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
//   Air->AddElement(N, 70*perCent);
//   Air->AddElement(O, 30*perCent);
		 
  //Copper
  G4Material* Cu = 
  new G4Material("Copper", z=29., a= 63.55*g/mole, density= 8.92*g/cm3);

//   //Aluminium for testing
//   G4Material* Cu = 
//   new G4Material("Aluminium", z=13., a= 26.98*g/mole, density= 2.7*g/cm3);
  
  //Silicon
  G4Material* Si = 
  new G4Material("Silicon", z=14., a= 28.09*g/mole, density= 2.33*g/cm3);

  // Vacuum
  density = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  a = 1.01*g/mole;
  z = 1;
  vacuum = new G4Material("Galactic", z, a,
			  density, kStateGas, temperature, pressure);

  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//--------- Sizes of the principal geometrical components (solids)  ---------
  
  // world volume
  hWorldLength   = 50.0*mm;

  // target
  hTargetLength  = 2.5*mm;
  hTargetDepth   = 15.0*mm;

  // CCD
  hTrackerLength = 20.0*mm;
  hTrackerDepth = 0.225*mm;

//--------- positions of the principal geometrical components --------

  // target
  xPosTarget = 0.0*cm;
  yPosTarget = 0.0*cm;
  zPosTarget = 0.0*cm;

  // tracker
  xPosTracker = 0.0*cm;
  yPosTracker = 2.0*cm;
  zPosTracker = 0.0*cm;

//--------- material names of the principal geometrical components --------
  
//   WorldMater   = Air;
  WorldMater   = vacuum;
  TargetMater  = Cu;
  TrackerMater = Si;
      
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  

  //------------------------------ 
  // World
  //------------------------------ 

  solidWorld= new G4Box("world",hWorldLength,hWorldLength,hWorldLength);
  logicWorld= new G4LogicalVolume(solidWorld,WorldMater,"World",0,0,0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // copy number

			 
  //------------------------------ 
  // Target
  //------------------------------
  
  G4ThreeVector positionTarget = 
    G4ThreeVector(xPosTarget,yPosTarget,zPosTarget);
   
  solidTarget = new G4Box("target",hTargetLength,hTargetLength,hTargetDepth);
  logicTarget = new G4LogicalVolume(solidTarget,TargetMater,"Target",0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  logicTarget,     // its logical volume
				  "Target",        // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // copy number 


  //------------------------------ 
  // Tracker
  //------------------------------
  
  G4ThreeVector positionTracker = 
    G4ThreeVector(xPosTracker,yPosTracker,zPosTracker);
  
  solidTracker = 
    new G4Box("tracker",hTrackerLength,hTrackerDepth,hTrackerLength);
  logicTracker = 
    new G4LogicalVolume(solidTracker,TrackerMater,"Tracker",0,0,0);  
  physiTracker = new G4PVPlacement(0,              // no rotation
				  positionTracker, // at (x,y,z)
				  logicTracker,    // its logical volume
				  "Tracker",       // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // copy number 




  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "eRosita/TrackerChamberSD";
  eRositaTrackerSD* aTrackerSD = new eRositaTrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  logicTracker->SetSensitiveDetector( aTrackerSD );


//--------- Visualization attributes -------------------------------

  // use this to make world volume invisible
  visWorld = new G4VisAttributes();
  visWorld->SetVisibility(false);
  logicWorld->SetVisAttributes(visWorld);

  // render target in redish color
  visTarget = new G4VisAttributes();
//   visTarget->SetColor(G4Color(1.0,0.3,0.3));  // redish
  visTarget->SetColor(G4Color(1.0,1.0,1.0));  // black
  logicTarget->SetVisAttributes(visTarget);

  // render tracker in blueish color
  visTracker = new G4VisAttributes();
//   visTracker->SetColor(G4Color(0.3,0.3,1.0));   // blueish
  visTracker->SetColor(G4Color(1.0,1.0,1.0));   // black
  logicTracker->SetVisAttributes(visTracker);


  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void eRositaDetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {TargetMater = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
             << materialName << G4endl;
     }             
}
*/ 

