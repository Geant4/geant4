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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProductionCuts.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
    // create commands for interactive definition of the detector
    fDetectorMessenger = new DetectorMessenger(this);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstruction::~DetectorConstruction()
{
    delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::Construct()

{
    DefineMaterials();
    return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstruction::DefineMaterials()
{ 
// Water is defined from NIST material database
    G4NistManager * man = G4NistManager::Instance();
    G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");

// Default materials in setup.
    fpWaterMaterial = H2O;

// needed variables

    G4double z, a, density;
    G4String name, symbol;
    G4int nComponents, nAtoms;


    a = 12.0107*g/mole;
    G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

    a = 1.00794*g/mole;
    G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z=1., a);

    a = 15.9994*g/mole;
    G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

    a = 14.0067*g/mole;
    G4Element* elN  = new G4Element(name="Nitrogen"  ,symbol="N" , z= 7., a);

// Definition of Tetrahydrofurane (THF)

    density=1.346*g/cm3;

    fpTHFMaterial = new G4Material("THF", density, nComponents=3);
    fpTHFMaterial->AddElement(elC, nAtoms=4);
    fpTHFMaterial->AddElement(elH, nAtoms=8);
    fpTHFMaterial->AddElement(elO, nAtoms=1);

// Definition of Nitrogen in nanodosimetry experiments

    density = 0.34e-6*g/cm3;

    fpN2Material = new G4Material("N2", density, nComponents=1);
    fpN2Material->AddElement(elN, nAtoms=2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
    G4double diameter;
    G4double highz;
    G4Material* targetMaterial;
    G4Material* worldMaterial;

    if (fGeomType == "dna"){
        // nanometric geometry
        fWorldSize  = 20.*nm;
        diameter = 2.3*nm;
        highz = 3.4*nm;
        targetMaterial = fpTHFMaterial;
        worldMaterial = fpWaterMaterial;
    } else {
        // macrometric geometry (experimental target in nanodosimetry)
        fWorldSize  = 2*cm;
        diameter = 1*cm;
        highz = 1*cm;;
        targetMaterial = fpN2Material;
        worldMaterial = fpN2Material;
    }


    fpSolidWorld = new G4Box("World", //its name
                             fWorldSize/2,fWorldSize/2,fWorldSize/2);  //its size

    fpLogicWorld = new G4LogicalVolume(fpSolidWorld,        //its solid
                                       worldMaterial,        //its material
                                       "World");                //its name

    fpPhysiWorld = new G4PVPlacement(0,                        //no rotation
                                     G4ThreeVector(),        //at (0,0,0)
                                     "World",                //its name
                                     fpLogicWorld,                //its logical volume
                                     0,                        //its mother  volume
                                     false,                        //no boolean operation
                                     0);                        //copy number

    G4Tubs* solidTarget = new G4Tubs("Target",0, diameter/2., highz/2., 0*degree, 360*degree);

    G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget,   //its solid
                                                       targetMaterial,  //its material
                                                       "Target");   //its name

    new G4PVPlacement(0,                             //no rotation
                      G4ThreeVector(),         //at (0,0,0)
                      "Target",                         //its name
                      logicTarget,                 //its logical volume
                      fpPhysiWorld,                 //its mother  volume
                      false,                         //no boolean operation
                      0);                             //copy number

// Visualization attributes
    G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    worldVisAtt->SetVisibility(true);
    fpLogicWorld->SetVisAttributes(worldVisAtt);

    G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    worldVisAtt1->SetVisibility(true);
    logicTarget->SetVisAttributes(worldVisAtt1);

    return fpPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void DetectorConstruction::SetGeometry(const G4String& name)
  {
    fGeomType = name;

    // tell RunManager about changes
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

  }
