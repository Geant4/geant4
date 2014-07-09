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

#include "Tst69DetectorConstruction.hh"

#include "Tst69DetectorMessenger.hh"

#include "G4SystemOfUnits.hh"
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
#include "G4Pow.hh"

Tst69DetectorConstruction::Tst69DetectorConstruction()
  :simpleBoxLog(0),selectedMaterial(0),theH(0),theLi(0),theC(0),theSi(0),theCu(0),thePb(0),theU(0),theTh(0),theMix(0)
{
  detectorMessenger = new Tst69DetectorMessenger(this);
  materialChoice = "mix";
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

  if(!theLi)
  {
    a = 7.*g/mole;
    density = 0.53*g/cm3;
    theLi = new G4Material(name="Lithium", z=3., a, density);
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

  if(!theMix) {
    // A fictitious mixed material with weight fractions such that the mean
    // free path for inelastic collisions with high-energy particles is roughly
    // the same for all components
    G4Pow *g4pow = G4Pow::GetInstance();
    G4Element* elH  = new G4Element("Hydrogen", "H", 1., 1.*g/mole);
    G4double weightH = 0.5*g4pow->A13(elH->GetA());
    G4Element* elHe = new G4Element("Helium", "He", 2., 4.*g/mole);
    G4double weightHe = 1.*g4pow->A13(elHe->GetA());
    G4Element* elC  = new G4Element("Carbon", "C", 6., 12.*g/mole);
    G4double weightC = 1.*g4pow->A13(elC->GetA());
    G4Element* elSi  = new G4Element("Silicon", "Si", 14., 28.*g/mole);
    G4double weightSi = 1.*g4pow->A13(elSi->GetA());
    G4Element* elCu  = new G4Element("Copper", "Cu", 29., 63.*g/mole);
    G4double weightCu = 1.*g4pow->A13(elCu->GetA());
    G4Element* elPb  = new G4Element("Lead", "Pb", 82., 208.*g/mole);
    G4double weightPb = 1.*g4pow->A13(elPb->GetA());
    G4Element* elU  = new G4Element("Uranium", "U", 92., 238.*g/mole);
    G4double weightU = 1.*g4pow->A13(elU->GetA());

    const G4double totalWeight = weightH + weightHe + weightC + weightSi + weightCu + weightPb + weightU;
    weightH /= totalWeight;
    weightHe /= totalWeight;
    weightC /= totalWeight;
    weightSi /= totalWeight;
    weightCu /= totalWeight;
    weightPb /= totalWeight;
    weightU /= totalWeight;

    theMix = new G4Material("mix", 4.*g/cm3, 7);
    theMix->AddElement(elH, weightH);
    theMix->AddElement(elHe, weightHe);
    theMix->AddElement(elC, weightC);  
    theMix->AddElement(elSi, weightSi);  
    theMix->AddElement(elCu, weightCu);  
    theMix->AddElement(elPb, weightPb);  
    theMix->AddElement(elU, weightU);  
  }

  if(materialChoice=="H")
    selectedMaterial = theH;
  else if(materialChoice=="Li")
    selectedMaterial = theLi;
  else if(materialChoice=="C")
    selectedMaterial = theC;
  else if(materialChoice=="Si")
    selectedMaterial = theSi;
  else if(materialChoice=="Cu")
    selectedMaterial = theCu;
  else if(materialChoice=="U")
    selectedMaterial = theU;
  else if(materialChoice=="Th")
    selectedMaterial = theTh;
  else if(materialChoice=="mix")
    selectedMaterial = theMix;
  else
    selectedMaterial = thePb;

  if(simpleBoxLog)
  { simpleBoxLog->SetMaterial(selectedMaterial); }
}

G4VPhysicalVolume* Tst69DetectorConstruction::Construct()
{
  SelectMaterialPointer();

  G4Box * mySimpleBox = new G4Box("SBox",20.*cm, 20.*cm, 20.*cm);
  simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                      selectedMaterial,"SLog",0,0,0);
  G4VPhysicalVolume* simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                        "SPhys",simpleBoxLog,0,false,0);

  return simpleBoxDetector;
}

