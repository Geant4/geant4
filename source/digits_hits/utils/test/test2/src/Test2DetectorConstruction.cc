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
// 

#include "Test2DetectorConstruction.hh"
#include "Test2GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


Test2DetectorConstruction::Test2DetectorConstruction()
  :G4VUserDetectorConstruction(),fbConstructed(false),
  GEOM(0),SDC(0),PSC(0),fSensitivityType(0)
{;}

Test2DetectorConstruction::~Test2DetectorConstruction()
{
  if(GEOM) delete GEOM;
  if(SDC) delete SDC;
  if(PSC) delete PSC;
}

G4VPhysicalVolume* Test2DetectorConstruction::Construct()
{
  if(!fbConstructed)
  { 
    fbConstructed = true;
    DefineMaterials();
    SetupGeometry();
    switch ( fSensitivityType ){
    case 0:
      break;
    case 1:
      SetupPSDetectors();
      break;
    case 2:
      SetupSDDetectors();
      break;
    default:
      break;
    }
  }
  return fWorldPhys;
}

void Test2DetectorConstruction::DefineMaterials()
{ 
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  

  G4int ncomponents, natoms;
  G4double fractionmass;

  //
  // define Elements
  //

  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //
  // define a material from elements.   case 1: chemical molecule
  //
 
  density = 1.000*g/cm3;
  fWaterMat = new G4Material(name="Water", density, ncomponents=2);
  fWaterMat->AddElement(H, natoms=2);
  fWaterMat->AddElement(O, natoms=1);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290*mg/cm3;
  fAirMat = new G4Material(name="Air"  , density, ncomponents=2);
  fAirMat->AddElement(N, fractionmass=0.7);
  fAirMat->AddElement(O, fractionmass=0.3);
}

void Test2DetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",2.*m,2.*m,2.*m);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,fAirMat,"World");
  fWorldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",
                        0,false,0);
  
  //                               
  // 3D nested phantom
  //
  // parameters
  phantomSize[0] = 1.*m;
  phantomSize[1] = 1.*m;
  phantomSize[2] = 1.*m;
  nSegment[0] = 10;
  nSegment[1] = 10;
  nSegment[2] = 10;

  GEOM = 
    new Test2GeometryConstruction(phantomSize,fWaterMat,nSegment);
  fWorldPhys = GEOM->ConstructGeometry(fWorldPhys);

}

#include "G4SDManager.hh"
#include "Test2PhantomSD.hh"
void Test2DetectorConstruction::SetupSDDetectors() {
  G4LogicalVolume* phantomLogical = GEOM->GetSensitiveLogical();

  SDC = new Test2SDConstruction("MassWorldSD",nSegment);
  SDC->SetupSensitivity(phantomLogical);
}

#include "Test2PSConstruction.hh"

void Test2DetectorConstruction::SetupPSDetectors() {
  G4LogicalVolume* phantomLogical = GEOM->GetSensitiveLogical();

  PSC = new Test2PSConstruction("MassWorldPS",nSegment);
  PSC->SetupSensitivity(phantomLogical);
}

