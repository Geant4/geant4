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
// $Id: DetectorConstruction.cc 73446 2013-08-27 11:32:59Z gcosmo $
//
//
// --------------------------------------------------------------
//      GEANT 4 - DetectorConstruction class
// --------------------------------------------------------------
//
// Author: Witold POKORSKI (Witold.Pokorski@cern.ch)
//
// --------------------------------------------------------------

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction() :  
  Iron(0), Copper(0), Tungsten(0), Lead(0), Uranium(0), PbWO4(0),
  Polystyrene(0), LiquidArgon(0), 
  theAbsorberMaterial(0),
  logicAbsorber(0), physiAbsorber(0) {}


DetectorConstruction::~DetectorConstruction() {}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //------------------- materials ------------------------

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, pressure, temperature, fractionmass;
  G4String name, symbol;
  G4int nel, natoms;

  //--- elements

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H2", z=1., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N2", z=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O2", z=8., a);

  a = 183.85*g/mole;
  G4Element* elW = new G4Element(name="Tungsten", symbol="W", z=74., a);

  a = 207.19*g/mole;
  G4Element* elPb = new G4Element(name="Lead", symbol="Pb", z=82., a);

  //--- simple materials

  // Iron has a  X0 = 1.7585 cm  and  lambda_I = 16.760 cm.   
  density = 7.87*g/cm3;
  a = 55.85*g/mole;
  Iron = new G4Material(name="Iron", z=26., a, density);

  // Copper has a  X0 = 1.4353 cm  and  lambda_I = 15.056 cm.   
  density = 8.96*g/cm3;
  a = 63.55*g/mole;
  Copper = new G4Material(name="Copper", z=29., a, density);

  // Tungsten has a  X0 = 0.35 cm  and  lambda_I = 9.5855 cm. 
  density = 19.30*g/cm3;
  a = 183.85*g/mole;
  Tungsten = new G4Material(name="Tungsten", z=74., a, density);

  // Lead has a  X0 = 0.56120 cm  and  lambda_I = 17.092 cm.  
  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  Lead = new G4Material(name="Lead", z=82., a, density);

  // Uranium has a  X0 = 0.31662 cm  and  lambda_I = 10.501 cm.  
  density =  18.95*g/cm3;
  a = 238.03*g/mole;
  Uranium = new G4Material(name="Uranium", z=92., a, density);

  // Liquid Argon has a  X0 = 10.971 cm  and  lambda_I = 65.769 cm.  
  density = 1.4*g/cm3;
  a = 39.95*g/mole;
  LiquidArgon = new G4Material(name="LiquidArgon", z=18., a, density);

  //--- mixtures

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, 0.7);
  Air->AddElement(elO, 0.3);

  // 4-May-2006 : We rename "Vacuum" as "G4vacuum" to avoid
  //              problems with Flugg.
  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* G4vacuum = new G4Material(name="G4vacuum", density, nel=1,
                                        kStateGas, temperature, pressure);
  G4vacuum->AddMaterial(Air, fractionmass=1.);

  // Plastic scintillator tiles (used both in CMS hadron calorimeter
  // and ATLAS hadron barrel calorimeter): 
  //     X0 = 42.4 cm  and  lambda_I = 79.360 cm.  
  density = 1.032*g/cm3;
  Polystyrene = new G4Material(name="Polystyrene", density, nel=2);
  Polystyrene->AddElement(elC, natoms=19);
  Polystyrene->AddElement(elH, natoms=21);

   // PbWO4 CMS crystals. It has a  X0 = 0.89 cm  and  lambda_I = 22.4 cm. 
  density = 8.28*g/cm3;
  PbWO4 = new G4Material(name="PbWO4", density, nel=3);
  PbWO4->AddElement(elPb, natoms=1);
  PbWO4->AddElement(elW,  natoms=1);
  PbWO4->AddElement(elO,  natoms=4);

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

  experimentalHall_log = new G4LogicalVolume(experimentalHall_box, // solid 
                                             G4vacuum,             // material
                                             "expHall_log",        // name
                                             0,                    // field manager
                                             0,                    // sensitive detector
                                             0);                   // user limits

  experimentalHall_phys = new G4PVPlacement(0,                     // rotation
                                            G4ThreeVector(),       // translation
                                            "expHall",             // name
                                            experimentalHall_log,  // logical volume
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

  logicAbsorber = new G4LogicalVolume(solidAbsorber,       // solid 
                                      theAbsorberMaterial, // material
                                      "logicAbsorber",     // name
                                      0,                   // field manager
                                      0,                   // sensitive detector
                                      0);                  // user limits

  physiAbsorber = new G4PVPlacement(0,                     // rotation
                                    G4ThreeVector(),       // translation
                                    "physiAbsorber",       // its name
                                    logicAbsorber,         // logical volume
                                    experimentalHall_phys, // mother physical volume
                                    false,                 // boolean operation
                                    100);                  // copy number

  // --- Set default values    ***LOOKHERE***
  theAbsorberMaterial = Iron;
  //theAbsorberMaterial = Copper;
  //theAbsorberMaterial = Tungsten;
  //theAbsorberMaterial = Lead;
  //theAbsorberMaterial = Uranium;
  //theAbsorberMaterial = PbWO4;
  //theAbsorberMaterial = Polystyrene;
  //theAbsorberMaterial = LiquidArgon;
  
  logicAbsorber->SetMaterial( theAbsorberMaterial );
  
  PrintParameters();

  return experimentalHall_phys;
}


void DetectorConstruction::PrintParameters()
{
  G4cout << G4endl << G4endl
         << " ------  DetectorConstruction::PrintParameters() ------ " << G4endl
         << " Absorber Material = ";
  if ( theAbsorberMaterial ) {
    G4cout << theAbsorberMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << " -------------------------------------------------------- "
         << G4endl << G4endl;
}
