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
// $Id: Tst69DetectorConstruction.cc,v 1.1 2008-11-10 11:06:21 kaitanie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst69DetectorConstruction.hh"

#include "Tst69DetectorMessenger.hh"

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

Tst69DetectorConstruction::Tst69DetectorConstruction()
  :simpleBoxLog(0),selectedMaterial(0),theH(0),theC(0),theSi(0),theCu(0),thePb(0),theU(0),theTh(0)
{
  detectorMessenger = new Tst69DetectorMessenger(this);
  materialChoice = "Pb";
}

Tst69DetectorConstruction::~Tst69DetectorConstruction()
{
  delete detectorMessenger;
}

void Tst69DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << G4endl;
}

void Tst69DetectorConstruction::SelectMaterialPointer()
{
//--------- Material definition ---------

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  if(!theH)
  {
    a = 1.007940*g/mole;
    G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
    density = 1.29e-03*g/cm3;
    theH = new G4Material(name="H", density, nel=1);
    theH->AddElement(elH, 1.0);
  }

  if(!theC)
  {
    a = 12.011*g/mole;
    density = 2.1*g/cm3;
    theC = new G4Material(name="Carbon", z=6., a, density);
  }

  if(!theSi)
  {
    a = 28.0855*g/mole;
    density = 2.33*g/cm3;
    theSi = new G4Material(name="Silicon", z=14., a, density);
  }

  if(!theCu)
  {
    a = 63.546*g/mole;
    density = 8.96*g/cm3;
    theCu = new G4Material(name="Copper", z=29., a, density);
  }

  if(!thePb)
  {
    a = 207.19*g/mole;
    z = 82.;
    density = 11.35*g/cm3;
    thePb = new G4Material(name="Lead", z=82., a, density);
  }

  if(!theU)
  {
    a = 238.0289*g/mole;
    density = 18.95*g/cm3;
    theU = new G4Material(name="Uranium", z=92., a, density);
  }

  if(!theTh)
  {
    a = 232.0381*g/mole;
    density = 11.7*g/cm3;
    theTh = new G4Material(name="Thorium", z=90., a, density);
  }

  if(materialChoice=="H")
  { selectedMaterial = theH; }
  else if(materialChoice=="C")
  { selectedMaterial = theC; }
  else if(materialChoice=="Si")
  { selectedMaterial = theSi; }
  else if(materialChoice=="Cu")
  { selectedMaterial = theCu; }
  else if(materialChoice=="U")
    {selectedMaterial = theU; }
  else if(materialChoice=="Th")
    {selectedMaterial = theTh; }
  else
  { selectedMaterial = thePb; }

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

G4VPhysicalVolume* Tst69DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * mySimpleBox = new G4Box("SBox",200*cm, 200*cm, 200*cm);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

