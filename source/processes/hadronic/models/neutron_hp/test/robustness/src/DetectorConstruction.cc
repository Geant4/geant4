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
// $Id$
//

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

DetectorConstruction::DetectorConstruction()
:simpleBoxLog(0),selectedMaterial(0),Air(0),Al(0),Pb(0),U(0)
{
  detectorMessenger = new DetectorMessenger(this);
  materialChoice = "Pb";
}

DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
}

void DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << G4endl;
}

void DetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;


    G4cout << "Please enter Z and A of the material"<<G4endl;
    G4cin >> z >> a;
    a*=g/mole;
    density = 5*g/cm3;
    selectedMaterial = new G4Material(name="Current", z, a, density);

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * mySimpleBox = new G4Box("SBox",200*cm, 200*cm, 200*cm);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

