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

#include "Tst1DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4SDParticleFilter.hh"

#include "G4ios.hh"

Tst1DetectorConstruction::Tst1DetectorConstruction()
:constructed(false)
{;}

Tst1DetectorConstruction::~Tst1DetectorConstruction()
{;}

G4VPhysicalVolume* Tst1DetectorConstruction::Construct()
{
  if(!constructed)
  { 
    constructed = true;
    DefineMaterials();
    SetupGeometry();
    SetupDetectors();
  }
  return worldPhys;
}

void Tst1DetectorConstruction::DefineMaterials()
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
  water = new G4Material(name="Water", density, ncomponents=2);
  water->AddElement(H, natoms=2);
  water->AddElement(O, natoms=1);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //

  density = 1.290*mg/cm3;
  air = new G4Material(name="Air"  , density, ncomponents=2);
  air->AddElement(N, fractionmass=0.7);
  air->AddElement(O, fractionmass=0.3);
}

void Tst1DetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",2.*m,2.*m,2.*m);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,air,"World");
  worldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",
                        0,false,0);
  
  //                               
  // Phantom
  //  
  G4VSolid* phantomSolid = new G4Box("Calor",1.*m,1.*m,1.*m);
  G4LogicalVolume* phantomLogical = new G4LogicalVolume(phantomSolid,water,"Phantom");
  phantomPhys = new G4PVPlacement(0,G4ThreeVector(),phantomLogical,"Phantom",
                         worldLogical,false,0);
  //                                        
  // Visualization attributes
  //
//  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  phantomLogical->SetVisAttributes(simpleBoxVisAtt);
}

void Tst1DetectorConstruction::SetupDetectors()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  G4LogicalVolume* phantomLogical = phantomPhys->GetLogicalVolume();

  G4String filterName, particleName;
  G4SDParticleFilter* gammaFilter = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  G4SDParticleFilter* electronFilter = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
  G4SDParticleFilter* positronFilter = new G4SDParticleFilter(filterName="positronFilter",particleName="e+");

  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("MassWorld");

  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("eDep");
  det->RegisterPrimitive(primitive);
  primitive = new G4PSTrackLength("trackLengthGamma");
  primitive->SetFilter(gammaFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSNofStep("nStepGamma");
  primitive->SetFilter(gammaFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSTrackLength("trackLengthElec");
  primitive->SetFilter(electronFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSNofStep("nStepElec");
  primitive->SetFilter(electronFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSTrackLength("trackLengthPosi");
  primitive->SetFilter(positronFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSNofStep("nStepPosi");
  primitive->SetFilter(positronFilter);
  det->RegisterPrimitive(primitive);

  G4SDManager::GetSDMpointer()->AddNewDetector(det);
  phantomLogical->SetSensitiveDetector(det);
  G4SDManager::GetSDMpointer()->SetVerboseLevel(0);
}

