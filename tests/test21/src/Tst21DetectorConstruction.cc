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
// $Id: Tst21DetectorConstruction.cc,v 1.4 2006-07-04 09:27:20 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst21DetectorConstruction.hh"

#include "Tst21DetectorMessenger.hh"

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

Tst21DetectorConstruction::Tst21DetectorConstruction()
:simpleBoxLog(0),selectedMaterial(0),Air(0),He(0),Be(0),C(0),Al(0),Pb(0)
{
  detectorMessenger = new Tst21DetectorMessenger(this);
  materialChoice = "Pb";
}

Tst21DetectorConstruction::~Tst21DetectorConstruction()
{
  delete detectorMessenger;
}

void Tst21DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << G4endl;
}

void Tst21DetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  if(!Air)
  {
    a = 14.01*g/mole;
    G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
    a = 16.00*g/mole;
    G4Element* elO = new G4Element(name="Oxigen", symbol="O", iz=8., a);
    density = 1.29e-03*g/cm3;
    Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);
  }

  if(!He)
  {
    a = 4.0027*g/mole;
    density = 2.*g/cm3; // Unnatural, just for the test
    He = new G4Material(name="Helium", z=2., a, density);
  }

  if(!Be)
  {
    a = 9.012*g/mole;
    density = 1.848*g/cm3;
    Be = new G4Material(name="Beryllium", z=4., a, density);
  }

  if(!C)
  {
    a = 12.*g/mole;
    density = 2.27*g/cm3;
    C = new G4Material(name="Carbon", z=6., a, density);
  }

  if(!Al)
  {
    a = 26.98*g/mole;
    density = 2.7*g/cm3;
    Al = new G4Material(name="Aluminium", z=13., a, density);
  }

  if(!Pb)
  {
    a = 207.19*g/mole;
    density = 11.35*g/cm3;
    Pb = new G4Material(name="Lead", z=82., a, density);
  }

  if    (materialChoice=="Air") selectedMaterial = Air;
  else if(materialChoice=="He") selectedMaterial = He;
  else if(materialChoice=="Be") selectedMaterial = Be;
  else if(materialChoice=="C")  selectedMaterial = C;
  else if(materialChoice=="Al") selectedMaterial = Al;
  else                          selectedMaterial = Pb;

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

G4VPhysicalVolume* Tst21DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * mySimpleBox = new G4Box("SBox",200*cm, 200*cm, 200*cm);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

