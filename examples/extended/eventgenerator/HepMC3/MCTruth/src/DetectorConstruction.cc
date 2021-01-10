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
/// \file eventgenerator/HepMC/MCTruth/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - DetectorConstruction class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
//
// --------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction() : 
  G4VUserDetectorConstruction(), 
  fAbsorberMaterial(0) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //------------------- materials ------------------------

  //--- simple materials

  G4NistManager* nistManager = G4NistManager::Instance();

  // Iron has a  X0 = 1.7585 cm  and  lambda_I = 16.760 cm.   
  G4Material* iron = nistManager->FindOrBuildMaterial("G4_Fe");

  // Copper has a  X0 = 1.4353 cm  and  lambda_I = 15.056 cm.   
  G4Material* copper = nistManager->FindOrBuildMaterial("G4_Cu");
  
  // Tungsten has a  X0 = 0.35 cm  and  lambda_I = 9.5855 cm. 
  G4Material* tungsten = nistManager->FindOrBuildMaterial("G4_W");

  // Lead has a  X0 = 0.56120 cm  and  lambda_I = 17.092 cm.  
  G4Material* lead = nistManager->FindOrBuildMaterial("G4_Pb");

  // Uranium has a  X0 = 0.31662 cm  and  lambda_I = 10.501 cm.  
  G4Material* uranium = nistManager->FindOrBuildMaterial("G4_U");

  // Liquid Argon has a  X0 = 10.971 cm  and  lambda_I = 65.769 cm.  
  G4double a, z, density;
  density = 1.4*g/cm3;
  a = 39.95*g/mole;
  G4Material* liquidArgon = new G4Material("LiquidArgon", z=18., a, density);

  //--- mixtures

  G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
  
  // 4-May-2006 : We rename "Vacuum" as "G4vacuum" to avoid
  //              problems with Flugg.
  G4double pressure, temperature, fractionMass;
  G4int nel;
  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* vacuum = new G4Material("G4vacuum", density, nel=1,
                                      kStateGas, temperature, pressure);
  vacuum->AddMaterial(air, fractionMass=1.);

  // Plastic scintillator tiles (used both in CMS hadron calorimeter
  // and ATLAS hadron barrel calorimeter): 
  //     X0 = 42.4 cm  and  lambda_I = 79.360 cm.  
  G4int natoms;
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  density = 1.032*g/cm3;
  G4Material* polystyrene = new G4Material("Polystyrene", density, nel=2);
  polystyrene->AddElement(elC, natoms=19);
  polystyrene->AddElement(elH, natoms=21);

   // PbWO4 CMS crystals. It has a  X0 = 0.89 cm  and  lambda_I = 22.4 cm. 
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elW = nistManager->FindOrBuildElement("W");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  density = 8.28*g/cm3;
  G4Material* pbWO4 = new G4Material("PbWO4", density, nel=3);
  pbWO4->AddElement(elPb, natoms=1);
  pbWO4->AddElement(elW,  natoms=1);
  pbWO4->AddElement(elO,  natoms=4);

  //------------------- volumes --------------------------

  // --- experimental hall (world volume)
  //     beam line along z axis

  //***LOOKHERE***
  const G4double sizeExpHall =  4.0*m;     // For normal calorimeter
  //const G4double sizeExpHall = 10.0*m;     // For Scintillator calorimeter

  G4double expHall_x = sizeExpHall / 2.0;  // half dimension along x 
  G4double expHall_y = sizeExpHall / 2.0;  // half dimension along y
  G4double expHall_z = sizeExpHall / 2.0;  // half dimension along z

  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* experimentalHallLog 
    = new G4LogicalVolume(experimentalHall_box, // solid 
                          vacuum,               // material
                          "expHall_log",        // name
                          0,                    // field manager
                          0,                    // sensitive detector
                          0);                   // user limits

  G4VPhysicalVolume* experimentalHallPhys 
    = new G4PVPlacement(0,                     // rotation
                        G4ThreeVector(),       // translation
                        "expHall",             // name
                        experimentalHallLog,   // logical volume
                        0,                     // mother physical volume
                        false,                 // boolean operation
                        0);                    // copy number
  
  // --- Detector

  //***LOOKHERE***
  const G4double sizeCalo = 2.0*m;         // For normal calorimeter
  //const G4double sizeCalo = 8.0*m;         // For Scintillator calorimeter

  G4double xAbsorber = sizeCalo / 2.0;  // half dimension along x 
  G4double yAbsorber = sizeCalo / 2.0;  // half dimension along y
  G4double zAbsorber = sizeCalo / 2.0;  // half dimension along z

  G4Box* solidAbsorber = new G4Box("solidAbsorber", xAbsorber, yAbsorber, zAbsorber);

  G4LogicalVolume* logicAbsorber 
    = new G4LogicalVolume(solidAbsorber,       // solid 
                          fAbsorberMaterial,   // material
                          "logicAbsorber",     // name
                          0,                   // field manager
                          0,                   // sensitive detector
                          0);                  // user limits

  new G4PVPlacement(0,                     // rotation
                    G4ThreeVector(),       // translation
                    "physiAbsorber",       // its name
                    logicAbsorber,         // logical volume
                    experimentalHallPhys,  // mother physical volume
                    false,                 // boolean operation
                    100);                  // copy number

  // --- Check if all materials were built
  if ( (! iron) || (! copper) || (! tungsten) || (! lead) || (! uranium) ||
       (! pbWO4) || (! polystyrene) || (! liquidArgon) ) {
    G4cerr << "Failure in building materials." << G4endl;
  }

  // --- Set default values    ***LOOKHERE***
  fAbsorberMaterial = iron;
  //fAbsorberMaterial = copper;
  //fAbsorberMaterial = tungsten;
  //fAbsorberMaterial = lead;
  //fAbsorberMaterial = uranium;
  //fAbsorberMaterial = pbWO4;
  //fAbsorberMaterial = polystyrene;
  //fAbsorberMaterial = liquidArgon;
  
  logicAbsorber->SetMaterial( fAbsorberMaterial );
  
  PrintParameters();

  return experimentalHallPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::PrintParameters()
{
  G4cout << G4endl << G4endl
         << " ------  DetectorConstruction::PrintParameters() ------ " << G4endl
         << " Absorber Material = ";
  if ( fAbsorberMaterial ) {
    G4cout << fAbsorberMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << " -------------------------------------------------------- "
         << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
