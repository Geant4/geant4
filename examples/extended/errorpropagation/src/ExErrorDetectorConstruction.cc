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

#include "ExErrorDetectorConstruction.hh"
#include "ExErrorDetectorMessenger.hh"
#include "ExErrorMagneticField.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include "G4Colour.hh"

#include "G4ios.hh"

//-------------------------------------------------------------
ExErrorDetectorConstruction::ExErrorDetectorConstruction()
  : xBEAM(5.*cm), xCDET(20.*cm), xECAL(40.*cm), xSOLN(10.*cm), xHCAL(100.*cm), 
    xMUON(50.*cm), ndivECAL(40./10.), ndivHCAL(100./10.), yzLength(50.*cm), 
    xHalfWorldLength(xBEAM + xCDET + xECAL + xSOLN + xHCAL + xMUON)
{

  // create UserLimits
  userLimits = new G4UserLimits();

  fpMagField = new ExErrorMagneticField(G4ThreeVector(0.*kilogauss,0.*kilogauss,-1.*kilogauss));
  detectorMessenger = new ExErrorDetectorMessenger(this);

}


//-------------------------------------------------------------
ExErrorDetectorConstruction::~ExErrorDetectorConstruction()
{
  delete fpMagField;
  delete detectorMessenger;             
}


//-------------------------------------------------------------
G4VPhysicalVolume* ExErrorDetectorConstruction::Construct()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  
  //Vacuum
  /*  a = 1.*g/mole;
  density = 1.E-9*g/cm3;
  G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a, density);
  */

  //Air
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
  density = 1.205*mg/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
  //Al
  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Al", z=13., a, density);
  //Fe
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Fe", z=26., a, density);
  //Cu
  a = 63.54*g/mole;
  density = 8.96*g/cm3;
  G4Material* Cu = new G4Material(name="Cu", z=29., a, density);

    
  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  //--------- Sizes of the principal geometrical components (solids)  --------- (half lengths)
  //double xBEAM  = 5.*2.*cm;
  //double xCDET  = 90.*cm;
  //double xECAL  = 40.*cm;
  //double xSOLN  = 10.*cm;
  //double xHCAL  = 100.*cm;
  //double xMUON  = 50.*cm;
  //double ndivECAL  = 10;
  //double ndivHCAL  = 10;
  //double yzLength = 100.*cm;
    
  //  double xWorldLength= xBEAM + xCDET + xECAL + xSOLN + xHCAL + xMUON;
  
   
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 
  //-  G4double HalfWorldLength = xWorldLength;
  G4cout << " HalfWorldLength " << xHalfWorldLength << G4endl;
  
  G4Box* solidWorld= new G4Box("world",xHalfWorldLength,yzLength,yzLength);
  G4LogicalVolume* logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
				 "World",         // its name
                                 logicWorld,      // its logical volume
                                 0,               // its mother  volume
                                 false,           // no boolean operations
                                 0);              // no field specific to volum
				 
  //------------------------------ 
  // BEAM
  //------------------------------   
  G4Box* solidBEAM = new G4Box("BEAM",xBEAM,yzLength,yzLength);
  G4LogicalVolume* logicBEAM = new G4LogicalVolume(solidBEAM,Air,"BEAM",0,0,0);
  G4ThreeVector positionBEAM = G4ThreeVector(0.,0.,0.);
  //G4VPhysicalVolume* physiBEAM = 
  new G4PVPlacement(0,               // no rotation
				  positionBEAM,  // at (x,y,z)
				  "BEAM",        // its name
				  logicBEAM,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //------------------------------ 
  // CDET (Central DETector)
  //------------------------------   
  G4ThreeVector positionCdet = G4ThreeVector(xBEAM + xCDET/2.,0.,0.);
  G4Box* solidCDET = new G4Box("CDET",xCDET/2.,yzLength,yzLength);
  G4LogicalVolume* logicCDET = new G4LogicalVolume(solidCDET,Air,"Cdet",0,0,0);
  //  G4VPhysicalVolume* physiCDET = 
  new G4PVPlacement(0,               // no rotation
				  positionCdet,  // at (x,y,z)
				  "CDET",        // its name
				  logicCDET,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //------------------------------ 
  // ECAL
  //------------------------------   
  G4ThreeVector positionECAL = G4ThreeVector(xBEAM + xCDET + xECAL/2., 0., 0.);
  G4Box* solidECAL = new G4Box("ECAL",xECAL/2.,yzLength,yzLength);
  G4LogicalVolume* logicECAL = new G4LogicalVolume(solidECAL,Cu,"ECAL",0,0,0);
  G4VPhysicalVolume* physiECAL = new G4PVPlacement(0,               // no rotation
				  positionECAL,  // at (x,y,z)
				  "ECAL",        // its name
				  logicECAL,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 
  //--------- Divide it 
  G4Box* solidECALdiv = new G4Box("ECAL",xECAL/2./ndivECAL,yzLength,yzLength);
  G4LogicalVolume* logicECALdiv = new G4LogicalVolume(solidECALdiv,Cu,"ECALdiv",0,0,0);
  new G4PVReplica("DVEC", logicECALdiv, physiECAL,
		      kXAxis, G4int(ndivECAL), xECAL/ndivECAL);


  //------------------------------ 
  // SOLN
  //------------------------------   
  G4ThreeVector positionSOLN = G4ThreeVector(xBEAM + xCDET + xECAL + xSOLN/2., 0., 0.);
  G4Box* solidSOLN = new G4Box("SOLN",xSOLN/2.,yzLength,yzLength);
  G4LogicalVolume* logicSOLN = new G4LogicalVolume(solidSOLN,Al,"SOLN",0,0,0);
  new G4PVPlacement(0,               // no rotation
				  positionSOLN,  // at (x,y,z)
				  "SOLN",        // its name
				  logicSOLN,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 

  //------------------------------ 
  // HCAL
  //------------------------------   
  G4ThreeVector positionHCAL = G4ThreeVector(xBEAM + xCDET + xECAL + xSOLN + xHCAL/2., 0., 0.);
  G4Box* solidHCAL = new G4Box("HCAL",xHCAL/2.,yzLength,yzLength);
  G4LogicalVolume* logicHCAL = new G4LogicalVolume(solidHCAL,Fe,"HCAL",0,0,0);
  G4VPhysicalVolume* physiHCAL = new G4PVPlacement(0,               // no rotation
				  positionHCAL,  // at (x,y,z)
				  "HCAL",        // its name
				  logicHCAL,     // its logical volume
				  physiWorld,      // its mother  volume
				  false,           // no boolean operations
				  0);              // no particular field 
  //--------- Divide it 
  G4Box* solidHCALdiv = new G4Box("HCAL",xHCAL/2./ndivHCAL,yzLength,yzLength);
  G4LogicalVolume* logicHCALdiv = new G4LogicalVolume(solidHCALdiv,Fe,"HCALdiv",0,0,0);
  new G4PVReplica("DVEH", logicHCALdiv, physiHCAL,
		  kXAxis, G4int(ndivHCAL), xHCAL/ndivHCAL);
  
  //------------------------------ 
  // MUON
  //------------------------------   
  G4ThreeVector positionMUON = G4ThreeVector(xBEAM + xCDET + xECAL + xSOLN + xHCAL + xMUON/2., 0., 0.);
  G4Box* solidMUON = new G4Box("MUON",xMUON/2.,yzLength,yzLength);
  G4LogicalVolume* logicMUON = new G4LogicalVolume(solidMUON,Air,"MUON",0,0,0);
  new G4PVPlacement(0,               // no rotation
		    positionMUON,  // at (x,y,z)
		    "MUON",        // its name
		    logicMUON,     // its logical volume
		    physiWorld,      // its mother  volume
		    false,           // no boolean operations
		    0);              // no particular field 


  G4VisAttributes* worldVisAtt = new G4VisAttributes(0);
  logicWorld->SetVisAttributes( worldVisAtt);
  return physiWorld;
}

 
//-------------------------------------------------------------
void ExErrorDetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetFieldValue(fieldValue);
}

