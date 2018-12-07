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
/// \file ExErrorDetectorConstruction.cc
/// \brief Implementation of the ExErrorDetectorConstruction class
//

#include "ExErrorDetectorConstruction.hh"
#include "ExErrorDetectorMessenger.hh"
#include "ExErrorMagneticField.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExErrorDetectorConstruction::ExErrorDetectorConstruction()
  : G4VUserDetectorConstruction(),
    fXBEAM(5.*cm), fXCDET(20.*cm), fXECAL(40.*cm), fXSOLN(10.*cm), fXHCAL(100.*cm), 
    fXMUON(50.*cm), fNdivECAL(40./10.), fNdivHCAL(100./10.), fYZLength(50.*cm), 
    fXHalfWorldLength(fXBEAM + fXCDET + fXECAL + fXSOLN + fXHCAL + fXMUON),
    fUserLimits(0), fMagField(0), fDetectorMessenger(0) 
{

  // create UserLimits
  fUserLimits = new G4UserLimits();

  fMagField = new ExErrorMagneticField(
                    G4ThreeVector(0.*kilogauss,0.*kilogauss,-1.*kilogauss));
  fDetectorMessenger = new ExErrorDetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExErrorDetectorConstruction::~ExErrorDetectorConstruction()
{
  delete fMagField;
  delete fDetectorMessenger;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* ExErrorDetectorConstruction::Construct()
{
//--------- Material definition ---------

  //Vacuum
  /*  a = 1.*g/mole;
  density = 1.E-9*g/cm3;
  G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a, density);
  */

  G4NistManager* nistMgr = G4NistManager::Instance();
  G4Material* air = nistMgr->FindOrBuildMaterial("G4_AIR");
  //Al
  G4Material* al = nistMgr->FindOrBuildMaterial("G4_Al");
  //Fe
  G4Material* fe = nistMgr->FindOrBuildMaterial("G4_Fe");
  //Cu
  G4Material* cu = nistMgr->FindOrBuildMaterial("G4_Cu");
    
  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  //--- Sizes of the principal geometrical components (solids)  --- (half lengths)
  //double fXBEAM  = 5.*2.*cm;
  //double fXCDET  = 90.*cm;
  //double fXECAL  = 40.*cm;
  //double fXSOLN  = 10.*cm;
  //double fXHCAL  = 100.*cm;
  //double fXMUON  = 50.*cm;
  //double fNdivECAL  = 10;
  //double fNdivHCAL  = 10;
  //double fYZLength = 100.*cm;
    
  //  double fXWorldLength= fXBEAM + fXCDET + fXECAL + fXSOLN + fXHCAL + fXMUON;
  
   
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
  
  //------------------------------ 
  // World
  //------------------------------ 
  //-  G4double HalfWorldLength = fXWorldLength;
  G4cout << " HalfWorldLength " << fXHalfWorldLength << G4endl;
  
  G4Box* solidWorld= new G4Box("world",fXHalfWorldLength,fYZLength,fYZLength);
  G4LogicalVolume* logicWorld= new G4LogicalVolume( solidWorld, air, "World", 0, 0, 0);
  //  Must place the World Physical volume unrotated at (0,0,0).
  // 
  G4VPhysicalVolume* physiWorld 
    = new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        "World",         // its name
                        logicWorld,      // its logical volume
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0);              // no field specific to volum
                                 
  //------------------------------ 
  // BEAM
  //------------------------------   
  G4Box* solidBEAM = new G4Box("BEAM",fXBEAM,fYZLength,fYZLength);
  G4LogicalVolume* logicBEAM = new G4LogicalVolume(solidBEAM,air,"BEAM",0,0,0);
  G4ThreeVector positionBEAM = G4ThreeVector(0.,0.,0.);
  //G4VPhysicalVolume* physiBEAM = 
  new G4PVPlacement(0,             // no rotation
                    positionBEAM,  // at (x,y,z)
                    "BEAM",        // its name
                    logicBEAM,     // its logical volume
                    physiWorld,    // its mother  volume
                    false,         // no boolean operations
                    0);            // no particular field 

  //------------------------------ 
  // CDET (Central DETector)
  //------------------------------   
  G4ThreeVector positionCdet = G4ThreeVector(fXBEAM + fXCDET/2.,0.,0.);
  G4Box* solidCDET = new G4Box("CDET",fXCDET/2.,fYZLength,fYZLength);
  G4LogicalVolume* logicCDET = new G4LogicalVolume(solidCDET,air,"Cdet",0,0,0);
  //  G4VPhysicalVolume* physiCDET = 
  new G4PVPlacement(0,             // no rotation
                    positionCdet,  // at (x,y,z)
                    "CDET",        // its name
                    logicCDET,     // its logical volume
                    physiWorld,    // its mother  volume
                    false,         // no boolean operations
                    0);            // no particular field 

  //------------------------------ 
  // ECAL
  //------------------------------   
  G4ThreeVector positionECAL = G4ThreeVector(fXBEAM + fXCDET + fXECAL/2., 0., 0.);
  G4Box* solidECAL = new G4Box("ECAL",fXECAL/2.,fYZLength,fYZLength);
  G4LogicalVolume* logicECAL = new G4LogicalVolume(solidECAL,cu,"ECAL",0,0,0);
  G4VPhysicalVolume* physiECAL 
    = new G4PVPlacement(0,             // no rotation
                        positionECAL,  // at (x,y,z)
                        "ECAL",        // its name
                        logicECAL,     // its logical volume
                        physiWorld,    // its mother  volume
                        false,         // no boolean operations
                        0);            // no particular field 
  //--------- Divide it 
  G4Box* solidECALdiv = new G4Box("ECAL",fXECAL/2./fNdivECAL,fYZLength,fYZLength);
  G4LogicalVolume* logicECALdiv = new G4LogicalVolume(solidECALdiv,cu,"ECALdiv",0,0,0);
  new G4PVReplica("DVEC", logicECALdiv, physiECAL,
                      kXAxis, G4int(fNdivECAL), fXECAL/fNdivECAL);


  //------------------------------ 
  // SOLN
  //------------------------------   
  G4ThreeVector positionSOLN 
    = G4ThreeVector(fXBEAM + fXCDET + fXECAL + fXSOLN/2., 0., 0.);
  G4Box* solidSOLN = new G4Box("SOLN",fXSOLN/2.,fYZLength,fYZLength);
  G4LogicalVolume* logicSOLN = new G4LogicalVolume(solidSOLN,al,"SOLN",0,0,0);
  new G4PVPlacement(0,             // no rotation
                    positionSOLN,  // at (x,y,z)
                    "SOLN",        // its name
                    logicSOLN,     // its logical volume
                    physiWorld,    // its mother  volume
                    false,         // no boolean operations
                    0);            // no particular field 

  //------------------------------ 
  // HCAL
  //------------------------------   
  G4ThreeVector positionHCAL = G4ThreeVector(fXBEAM + fXCDET + fXECAL + fXSOLN + fXHCAL/2., 0., 0.);
  G4Box* solidHCAL = new G4Box("HCAL",fXHCAL/2.,fYZLength,fYZLength);
  G4LogicalVolume* logicHCAL = new G4LogicalVolume(solidHCAL,fe,"HCAL",0,0,0);
  G4VPhysicalVolume* physiHCAL 
    = new G4PVPlacement(0,             // no rotation
                        positionHCAL,  // at (x,y,z)
                        "HCAL",        // its name
                        logicHCAL,     // its logical volume
                        physiWorld,    // its mother  volume
                        false,         // no boolean operations
                        0);            // no particular field 
  //--------- Divide it 
  G4Box* solidHCALdiv = new G4Box("HCAL",fXHCAL/2./fNdivHCAL,fYZLength,fYZLength);
  G4LogicalVolume* logicHCALdiv = new G4LogicalVolume(solidHCALdiv,fe,"HCALdiv",0,0,0);
  new G4PVReplica("DVEH", logicHCALdiv, physiHCAL,
                  kXAxis, G4int(fNdivHCAL), fXHCAL/fNdivHCAL);
  
  //------------------------------ 
  // MUON
  //------------------------------   
  G4ThreeVector positionMUON 
    = G4ThreeVector(fXBEAM + fXCDET + fXECAL + fXSOLN + fXHCAL + fXMUON/2., 0., 0.);
  G4Box* solidMUON = new G4Box("MUON",fXMUON/2.,fYZLength,fYZLength);
  G4LogicalVolume* logicMUON = new G4LogicalVolume(solidMUON,air,"MUON",0,0,0);
  new G4PVPlacement(0,             // no rotation
                    positionMUON,  // at (x,y,z)
                    "MUON",        // its name
                    logicMUON,     // its logical volume
                    physiWorld,    // its mother  volume
                    false,         // no boolean operations
                    0);            // no particular field 


  G4VisAttributes* worldVisAtt = new G4VisAttributes(0);
  logicWorld->SetVisAttributes( worldVisAtt);
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExErrorDetectorConstruction::SetMagField(G4double fieldValue)
{
  fMagField->SetFieldValue(fieldValue);
}

