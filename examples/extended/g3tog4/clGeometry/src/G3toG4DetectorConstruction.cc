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
/// \file g3tog4/clGeometry/src/G3toG4DetectorConstruction.cc
/// \brief Implementation of the G3toG4DetectorConstruction class
//
//
//
//--------------------------------------------------------------------------
// G3toG4DetectorConstruction. Most the work is Done in
// G4BuildGeom, which returns a G4LogicalVolume*, a pointer to the
// top-level logiical volume in the detector defined by the call List file
// inFile
//--------------------------------------------------------------------------

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G3toG4DetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G3toG4DetectorConstruction::G3toG4DetectorConstruction(G4String inFile)
  : G4VUserDetectorConstruction(),
    fInFile(inFile)
{
  G4cout << "Instantiated G3toG4DetectorConstruction using call list file \"" 
         << fInFile << "\"" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G3toG4DetectorConstruction::~G3toG4DetectorConstruction()
{
  // G4cout << "Deleted G3toG4DetectorConstruction..." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* G3toG4DetectorConstruction::Construct()
{
  G4LogicalVolume* lv = G4BuildGeom(fInFile);
  //G4LogicalVolume* lv = SimpleConstruct();
  if ( lv ) {
    G4VPhysicalVolume* pv
      = new G4PVPlacement(0, G4ThreeVector(), lv, lv->GetName(), 0, false, 0);
    G4cout << "Top-level G3toG4 logical volume " << lv->GetName() << " "
           << *(lv -> GetVisAttributes()) << G4endl;
    return pv;       
  } 

  G4cerr << "creation of logical mother failed !!!" << G4endl;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* G3toG4DetectorConstruction::SimpleConstruct()
{
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density, fractionmass; //z=mean number of protons;
  G4int ncomponents;          //iz=number of protons  in an isotope;
  // n=number of nucleons in an isotope;

  a = 14.01*g/mole;
  G4Element* eN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* eO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290*mg/cm3;
  G4Material* air = new G4Material(name="Air"  , density, ncomponents=2);
  air->AddElement(eN, fractionmass=0.7);
  air->AddElement(eO, fractionmass=0.3);
  G4VSolid* mother = new G4Box("TestMother",             //its name
                                100*cm, 100*cm, 100*cm); //its size

  G4VSolid* daughter = new G4Box("TestDaughter", 50*cm, 20*cm, 10*cm);

  G4LogicalVolume* logicMother = new G4LogicalVolume(mother,     //its solid
                                                     air,        //its material
                                                     "LTestMother");//its name

  G4LogicalVolume* logicDaughter = new G4LogicalVolume(daughter, //its solid
                                                       air,      //its material
                                                       "LTestDaughter");

  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicDaughter,
                    "PTestDaughter",
                    logicMother,
                    false,0);
  //
  // Visualization attributes
  //

  logicMother->SetVisAttributes (G4VisAttributes::GetInvisible());
  G4VisAttributes* daughterVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  daughterVisAtt->SetVisibility(true);
  logicDaughter->SetVisAttributes(daughterVisAtt);
  return logicMother;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
