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
// $Id: F07DetectorConstruction.cc,v 1.22 2010-01-22 11:57:03 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "F07DetectorConstruction.hh"
#include "F07DetectorMessenger.hh"
// #include "F07ChamberParameterisation.hh"
#include "F07MagneticField.hh"
#include "F07TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
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
 
F07DetectorConstruction::F07DetectorConstruction()
:solidWorld(0),  logicWorld(0),  physiWorld(0),
 solidTarget(0), logicTarget(0), physiTarget(0), 
 solidTracker(0),logicTracker(0),physiTracker(0), 
 fTargetMaterial(0), fChamberAbsorberMat(0), fChamberGasMat(0),
 stepLimit(0), fpMagField(0),
 fWorldLength(0.),  fTrackerLength(0.),  fTrackerRadius(0.),
 fNbOfChambers(0) , fTargetLength(0.),  fTargetMaxRadius(0.)
{
  fpMagField = new F07MagneticField();
  detectorMessenger = new F07DetectorMessenger(this);
  // solidChamber(0),logicChamber(0),physiChamber(0), 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
F07DetectorConstruction::~F07DetectorConstruction()
{
  delete fpMagField;
  delete stepLimit;
//  delete chamberParam;
  delete detectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* F07DetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel;
  G4bool checkPlacement= true;
    
  //Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);
   
  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  //Lead
  G4Material* Pb = 
  new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);
  G4Material* Fe= G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe"); 
    
  //Xenon gas
  G4Material* Xenon = 
  new G4Material("XenonGas", z=54., a=131.29*g/mole, density= 5.458*mg/cm3,
                 kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // Choice of material
  fTargetMaterial  = Pb;
  fChamberAbsorberMat = Fe;   // Absorber
  fChamberGasMat = Xenon;  // 

  G4cout << G4endl << " Chosen Materials: " << G4endl;
  G4cout << "  - Target          : " << fTargetMaterial->GetName() << G4endl;
  G4cout << "  - Chamber Absorber: " << fChamberAbsorberMat->GetName() << G4endl;
  G4cout << "  - Chamber Gas     : " << fChamberGasMat->GetName() << G4endl;


//--------- Sizes of the principal geometrical components (solids)  ---------
  
  fNbOfChambers = 20;      // Even is best: inner = gas, outermost = solid
  // ChamberWidth = 20*cm;
  // ChamberSpacing = 80*cm;

  fTrackerLength  = 5000.0 * mm;       // Full length of Tracker  
  fTrackerRadius  = 3000.0 * mm;       
  G4double fForwardSphereRadius = fTrackerRadius;

  G4double DeltaRadius=         fTrackerRadius/(fNbOfChambers+1);
  G4double TrackerLayerRadius=  DeltaRadius;
    
  fTargetLength    =  10.0 * mm;        // Full length           of Target
  fTargetMaxRadius =   1.0 * mm;        // Maximum Radial extent of Target
  fTargetMaxRadius = std::min( TrackerLayerRadius, fTargetMaxRadius ); 

  fWorldLength=     2.05 * std::max(fTrackerLength, fForwardSphereRadius);
  fWorldHalfWidth=  1.05 * std::max(fTrackerRadius, fForwardSphereRadius); 
   
  // G4double trackerSize = 0.5*fTrackerLength;   // Half length of the Tracker
      
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  solidWorld= new G4Box("world",fWorldHalfWidth,fWorldHalfWidth,0.5*fWorldLength);
  logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
                                 logicWorld,      // its logical volume
				 "World",         // its name
                                 0,               // its mother (logical) volume
                                 false,           // no boolean operations
                                 0);              // copy number
				 
  //------------------------------ 
  // Tracker 'Envelope' 
  //------------------------------
  G4double TrackerEnvInnerRadius= DeltaRadius;
  G4double TrackerEnvZhalfLength= (fTrackerLength+fForwardSphereRadius)/2.0;

  solidTracker = new G4Tubs("Tracker Chamber-Tubs", TrackerEnvInnerRadius, fTrackerRadius, 
			    TrackerEnvZhalfLength, 0.0,  360.0*degree); // Start & Delta Phi

  G4cout << "Created Tracker Envelope Cylinder: " << *solidTracker << G4endl;

  G4ThreeVector positionTracker = 
        G4ThreeVector(0.0, 
                      0.0, 
                      (fTrackerRadius-TrackerEnvZhalfLength) );
  
  logicTracker = new G4LogicalVolume(solidTracker , Air, "Tracker-LV"); // ,0,0,0);  
  physiTracker = new G4PVPlacement(0,              // no rotation
				  positionTracker, // at (x,y,z)
				  logicTracker,    // its logical volume				  
				  "TrackerEnvelope-PV",  // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  1);              // copy number 

  //---------------------------------------- 
  // Tracker Barrel - Cylindrical Layers
  //---------------------------------------- 
  // 
  // An example of Parameterised volumes
  // dummy values for G4Box -- modified by parameterised volume

  G4double ChamberInnerRadius, ChamberOuterRadius;
  G4double ChamberZhalfLength= 0.5 * fTrackerLength;

  unsigned int iChamber;
  for ( iChamber=0; iChamber<fNbOfChambers; iChamber++)
  {
     G4VSolid*          oneSolidChamber;
     G4Material*        ptrMaterial; 
     G4LogicalVolume*   oneLogicChamber; 
     G4VPhysicalVolume* physiChamber;
      
     ChamberInnerRadius= (iChamber+1)*DeltaRadius; 
     ChamberOuterRadius= (iChamber+2)*DeltaRadius; 
     oneSolidChamber = new G4Tubs("Tracker Chamber-Tubs", 
				  ChamberInnerRadius, 
				  ChamberOuterRadius,
				  ChamberZhalfLength, 
				  0.0,                // Start Phi
				  360.0*degree);      // Delta Phi
     // solidChambers.push(oneSolidChamber); 

     ptrMaterial= ((iChamber%2)==0) ? fChamberAbsorberMat : fChamberGasMat; // Mater->Absorber

     oneLogicChamber = new G4LogicalVolume(oneSolidChamber, ptrMaterial, "Chamber"); // ,0,0,0);
     // logicalChambers.push(oneLogicChamber); 

      G4ThreeVector positionChamber( 0.0, 0.0, 0.0 ); 
     physiChamber = new G4PVPlacement( 0,               // no rotation
				       positionChamber, // 
				       oneLogicChamber,    // its logical volume
				       "Chamber-PV",    // name
				       logicTracker,    // Mother logical volume
				       false,           // no boolean
				       iChamber,        // copy number
				       checkPlacement );   
     // physiChambers.push(physiChamber); 
     G4cout << " Tracker Chamber number " << iChamber << " placed " 
	    << " Radius = " << ChamberInnerRadius 
	    << " to " << ChamberOuterRadius << G4endl;
  }

  G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
         << "The chambers are made of layers of " << DeltaRadius/mm << " mm of " 
         << fChamberAbsorberMat->GetName() << " and "
         << fChamberGasMat->GetName() << "." << G4endl
         << "\n The length of the chambers is "
	 << fTrackerLength/mm << " mm"
    << " and their radius is " << fTrackerRadius/mm << "mm. "
    << G4endl;


  //------------------------------ 
  // Target
  //------------------------------
  G4ThreeVector positionTarget = G4ThreeVector(0,0,-(0.5*targetLength));

   //  The +z face of the tracker should have its center at (0.0, 0.0, 0.0)
  solidTarget = new G4Cons("target-Cone", 
			   0.0,        // Rmin 1
			   0.0,        // Rmax 1
			   0.0,        // Rmin 2
			   fTargetMaxRadius, // Rmax 2
			   targetHalfLength,  // Dz
			   0.0,        // Start Phi
			   360.0*degree); // Delta Phi
  // G4double targetSize  =  fTargetMaxRadius / std::sqrt(2.0); // Half width
  // new G4Box( "target-Box", targetSize, targetSize, targetHalfLength);
  logicTarget = new G4LogicalVolume(solidTarget,fTargetMaterial,"Target"); // ,0,0,0);
  physiTarget = new G4PVPlacement(0,               // no rotation
				  positionTarget,  // at (x,y,z)
				  logicTarget,     // its logical volume				  
				  "Target",        // its name
				  logicWorld,      // its mother  volume
				  false,           // no boolean operations
				  1,               // copy number 
				  true);           // check positioning

  G4cout << "Target is " << fTargetLength/cm << " cm of " 
         << fTargetMaterial->GetName() << G4endl;

	 
  //------------------------------------------------ 
  // Sensitive detectors
  //------------------------------------------------ 

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String trackerChamberSDname = "F07/TrackerChamberSD";
  F07TrackerSD* aTrackerSD = new F07TrackerSD( trackerChamberSDname );
  SDman->AddNewDetector( aTrackerSD );
  oneLogicChamber->SetSensitiveDetector( aTrackerSD );

//--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  logicWorld  ->SetVisAttributes(BoxVisAtt);  
  logicTarget ->SetVisAttributes(BoxVisAtt);
  logicTracker->SetVisAttributes(BoxVisAtt);
  
  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  oneLogicChamber->SetVisAttributes(ChamberVisAtt);
  
//--------- example of User Limits -------------------------------

  // below is an example of how to set tracking constraints in a given
  // logical volume(see also in F07PhysicsList how to setup the processes
  // G4StepLimiter or G4UserSpecialCuts).
    
  // Sets a max Step length in the tracker region, with G4StepLimiter
  //
  G4double maxStep = 0.5*ChamberWidth;
  stepLimit = new G4UserLimits(maxStep);
  logicTracker->SetUserLimits(stepLimit);
  
  // Set additional contraints on the track, with G4UserSpecialCuts
  //
  // G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  // logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
  //                                               minEkin));
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void F07DetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {fTargetMaterial = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The target is " << fTargetLength/cm << " cm of "
             << materialName << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F07DetectorConstruction::setChamberMaterial(G4String materialName)
{
  // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  
  if (pttoMaterial)
     {fChamberAbsorberMater = pttoMaterial;
      oneLogicChamber->SetMaterial(pttoMaterial); 
      G4cout << "\n----> The chambers will be created with " 
         << fDeltaRadius/mm << " mm of "
             << materialName << " as absorber. " << G4endl;
         G4cout<< " NOTE: existing chambers were NOT changed." << G4endl;
     }             
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void F07DetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetMagFieldValue(fieldValue);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F07DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
