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

#include "Tst11DetectorConstruction.hh"

#include "Tst11DetectorMessenger.hh"

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

Tst11DetectorConstruction::Tst11DetectorConstruction()
:simpleBoxLog(0),selectedMaterial(0),Air(0),Al(0),Pb(0),U(0)
{
  detectorMessenger = new Tst11DetectorMessenger(this);
  materialChoice = "Pb";
}

Tst11DetectorConstruction::~Tst11DetectorConstruction()
{
  delete detectorMessenger;
}

void Tst11DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << G4endl;
}

void Tst11DetectorConstruction::SelectMaterialPointer()
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

  if(!U)
  {
    a = 238.*g/mole;
    density = 21.*g/cm3;
    U  = new G4Material(name="Uranium", density, 1);
    G4Element * elU = new G4Element(name="Uranium", "U138", 1);
    G4Isotope * isoU = new G4Isotope("U238", 92, 238, a);
    elU->AddIsotope(isoU, 100.0);
    U->AddElement(elU, 1);
  }

  if(materialChoice=="Air")
  { selectedMaterial = Air; }
  else if(materialChoice=="Al")
  { selectedMaterial = Al; }
  else if(materialChoice=="Pb")
  { selectedMaterial = Pb; }
  else
  { selectedMaterial = U; }

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

G4VPhysicalVolume* Tst11DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * mySimpleBox = new G4Box("SBox",5*cm, 5*cm, 20*m);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

