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
#include "Tst28DetectorConstruction.hh"

#include "Tst28DetectorMessenger.hh"

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
#include "G4RunManager.hh"
#include "G4ios.hh"

Tst28DetectorConstruction::Tst28DetectorConstruction()
 : simpleBoxLog(0), simpleBoxDetector(0), theH(0), theSi(0), theCu(0), thePb(0),
   theU(0), theH2O(0), selectedMaterial(0)
{
  detectorMessenger = new Tst28DetectorMessenger(this);
  materialChoice = "Pb";
}

Tst28DetectorConstruction::~Tst28DetectorConstruction()
{
  delete detectorMessenger;
}

void Tst28DetectorConstruction::SelectMaterial(G4String val)
{
  materialChoice = val;
  SelectMaterialPointer();
  G4cout << "SimpleBox is now made of " << materialChoice << G4endl;
}

void Tst28DetectorConstruction::SelectMaterialPointer()
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
    a = 207.2*g/mole;
    density = 11.35*g/cm3;
    thePb = new G4Material(name="Lead", z=82., a, density);
  }

  if(!theU)
  {
    a = 238.0289*g/mole;
    density = 18.95*g/cm3;
    theU = new G4Material(name="Uranium", z=92., a, density);
  }

  if(!theH2O)
    {
    G4int ncomponents, Z, N;
    G4double abundance;

    // Oxygen
    G4Isotope* O16 = new G4Isotope(name="O16", Z=8, N=16, a = 15.994*g/mole);
    G4Isotope* O17 = new G4Isotope(name="O17", Z=8, N=17, a = 16.999*g/mole);
    G4Isotope* O18 = new G4Isotope(name="O18", Z=8, N=18, a = 17.999*g/mole);
    G4Element* elO  = new G4Element(name="Oxygen", symbol="O", ncomponents=3);
    elO->AddIsotope(O16, abundance= 99.762*perCent);
    elO->AddIsotope(O17, abundance= 0.038*perCent);
    elO->AddIsotope(O18, abundance= 0.200*perCent);

    // Hydrogen
    G4Isotope* H1 = new G4Isotope(name="H1", Z=1, N=1, a = 1.008*g/mole);
    G4Isotope* H2 = new G4Isotope(name="H2", Z=1, N=2, a = 2.014*g/mole);
    G4Element* elH  = new G4Element(name="Hydrogen", symbol="H", ncomponents=2);
    elH->AddIsotope(H1, abundance= 99.985*perCent);
    elH->AddIsotope(H2, abundance= 0.015*perCent);

    // H2O
    theH2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
    theH2O->AddElement(elH, 2);
    theH2O->AddElement(elO, 1);
    // overwrite computed meanExcitationEnergy with ICRU recommended value
    theH2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  }
  if(materialChoice=="H")
  { selectedMaterial = theH; }
  else if(materialChoice=="Si")
  { selectedMaterial = theSi; }
  else if(materialChoice=="Cu")
  { selectedMaterial = theCu; }
  else if(materialChoice=="Pb")
  { selectedMaterial = thePb; }
  else if(materialChoice=="H2O")
  { selectedMaterial = theH2O; }
    else
  { selectedMaterial = theU; }

  if(simpleBoxLog)
  {
    simpleBoxLog->SetMaterial(selectedMaterial);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

G4VPhysicalVolume* Tst28DetectorConstruction::Construct()
{
  SelectMaterialPointer();
  if(!simpleBoxLog) {
    G4Box * mySimpleBox = new G4Box("SBox",200*cm, 200*cm, 200*cm);
    simpleBoxLog = new G4LogicalVolume( mySimpleBox,
                                        selectedMaterial,"SLog",0,0,0);
    simpleBoxDetector = new G4PVPlacement(0,G4ThreeVector(),
                                          "SPhys",simpleBoxLog,0,false,0);
  }
  return simpleBoxDetector;
}

void Tst28DetectorConstruction::ConstructSDandField()
{}

