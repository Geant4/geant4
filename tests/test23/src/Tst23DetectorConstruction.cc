//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst23DetectorConstruction.cc,v 1.1 2001-12-14 14:53:41 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst23DetectorConstruction.hh"

#include "Tst23DetectorMessenger.hh"
#include "Tst23SensitiveDetector.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"



G4int Tst23DetectorConstruction::materialIndex =-1;
const G4String Tst23DetectorConstruction::materialName[N_MaterialType]
  = {  "Air", "Al", "Pb" };

Tst23DetectorConstruction::Tst23DetectorConstruction()
{
  for (size_t idx=0; idx<N_Regions; idx++) {
    simpleBoxLog[idx]=0;
    selectedMaterial[idx]=0;
  }
  detectorMessenger = new Tst23DetectorMessenger(this);
}

Tst23DetectorConstruction::~Tst23DetectorConstruction()
{
  delete detectorMessenger;
}

void Tst23DetectorConstruction::SelectMaterial(G4String val)
{
  for (size_t idx=0; idx<N_MaterialType; idx++) {
    if (val == materialName[idx]) {
      materialIndex = idx;
      SelectMaterialPointer();
      G4cout << "SimpleBox is now made of " << val << G4endl;
    }
  }
}

void Tst23DetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  static G4bool isInitialized = false;

  if (!isInitialized) {
    isInitialized = true;
    // Air
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
    density = 1.29e-03*g/cm3;
    // Material 
    materials[0][0] = new G4Material(materialName[0], density, nel=2);
    materials[0][0]->AddElement(elN, .7);
    materials[0][0]->AddElement(elO, .3);
    materials[0][1] = new G4Material(materialName[0]+"'", density, nel=2);
    materials[0][1]->AddElement(elN, .7);
    materials[0][1]->AddElement(elO, .3);

    // Al
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    materials[1][0] = new G4Material(materialName[1], z=13., a, density);
    materials[1][1] = new G4Material(materialName[1]+"'", z=13., a, density);
  
    // Pb
    a = 207.19*g/mole;
    density = 11.35*g/cm3;
    materials[2][0] = new G4Material(materialName[2], z=82., a, density);
    materials[2][1] = new G4Material(materialName[2]+"'", z=82., a, density);
  }

  selectedMaterial[0] = materials[materialIndex][0];
  selectedMaterial[1] = materials[materialIndex][1];

  if(simpleBoxLog[0]) { 
    simpleBoxLog[0]->SetMaterial(selectedMaterial[0]); 
    simpleBoxLog[1]->SetMaterial(selectedMaterial[1]); 
  }
}

G4VPhysicalVolume* Tst23DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * worldBox = new G4Box("WORLD",200*cm, 200*cm, 400*cm);
  G4LogicalVolume* worldLog = new G4LogicalVolume( worldBox,
                                      materials[0][0],"WLog",0,0,0);
  G4VPhysicalVolume* world = new G4PVPlacement(0,G4ThreeVector(),
                                        "World",worldLog,0,false,0);

  G4Box * mySimpleBox = new G4Box("SBox",200*cm, 200*cm, 200*cm);
  simpleBoxLog[0] = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial[0],"SLog0",0,0,0);
  simpleBoxLog[1] = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial[1],"SLog1",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector0 
    = new G4PVPlacement(0,
			G4ThreeVector(0.,0.,200*cm),
			"SPhys",
			simpleBoxLog[0],
			world,
			false,0);
  G4VPhysicalVolume* simpleBoxDetector1 
    = new G4PVPlacement(0,
			G4ThreeVector(0.,0.,-200*cm),
			"SPhys",
			simpleBoxLog[1],
			world,
			false,1);

  //--------- Sensitive detector -------------------------------------
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname = "Tst23";
  Tst23SensitiveDetector* pSD = new Tst23SensitiveDetector( SDname );
  SDman->AddNewDetector(  pSD );
   
  simpleBoxLog[0]->SetSensitiveDetector(pSD);
  simpleBoxLog[1]->SetSensitiveDetector(pSD);


  return world; 

} 

const G4String& Tst23DetectorConstruction::GetMaterialName(G4int idx)
{
  static G4String name;
  name = materialName[materialIndex];
  if (idx==1)   name += "'";

  return name;
}

