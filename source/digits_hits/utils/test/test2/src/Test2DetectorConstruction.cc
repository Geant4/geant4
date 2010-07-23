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
// $Id: Test2DetectorConstruction.cc,v 1.1 2010-07-23 06:15:41 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Test2DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

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

#include "Test2PhantomSD.hh"

Test2DetectorConstruction::Test2DetectorConstruction()
:constructed(false)
{;}

Test2DetectorConstruction::~Test2DetectorConstruction()
{;}

G4VPhysicalVolume* Test2DetectorConstruction::Construct()
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

void Test2DetectorConstruction::SetupGeometry()
{
  //     
  // World
  //
  G4VSolid* worldSolid = new G4Box("World",2.*m,2.*m,2.*m);
  G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,air,"World");
  worldPhys = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"World",
                        0,false,0);
  
  //                               
  // 3D nested phantom
  //
  // parameters
  G4double phantomSize[3] = {1.*m, 1.*m, 1.*m};
  G4int nSegment[3] = {10, 10, 10};

  //  phantom box
  G4VSolid* phantomSolid = new G4Box("PhantomBox", phantomSize[0], phantomSize[1], phantomSize[2]);
  G4LogicalVolume* phantomLogical = new G4LogicalVolume(phantomSolid, water, "Phantom");
  phantomPhys = new G4PVPlacement(0,G4ThreeVector(), phantomLogical, "Phantom",
                         worldLogical, false, 0);

  G4String layerName[3] = {"layerX", "layerY", "layerZ"};
  G4VSolid * layerSolid[3];
  G4LogicalVolume * layerLogical[3];

  //
  // replication along X
  //
  layerSolid[0] = new G4Box(layerName[0],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1],
			    phantomSize[2]);
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], water, layerName[0]);
  new G4PVReplica(layerName[0], layerLogical[0], phantomLogical, kXAxis,
		  nSegment[0], phantomSize[0]/nSegment[0]*2.);

  // replication along Y
  layerSolid[1] = new G4Box(layerName[1],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1]/nSegment[1],
			    phantomSize[2]);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], water, layerName[1]);
  new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], kYAxis,
		  nSegment[1], phantomSize[1]/nSegment[1]*2.);

  // replication along Z
  layerSolid[2] = new G4Box(layerName[2],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1]/nSegment[1],
			    phantomSize[2]/nSegment[2]);
  layerLogical[2] = new G4LogicalVolume(layerSolid[2], water, layerName[2]);
  new G4PVReplica(layerName[2], layerLogical[2], layerLogical[1], kZAxis,
		  nSegment[2], phantomSize[2]/nSegment[2]*2.);

  //
  // sensitive detectors
  //
  G4SDManager * sdMngr = G4SDManager::GetSDMpointer();

  G4String phantomSDName = "Test2/Phantom";
  Test2PhantomSD * phantomSD = new Test2PhantomSD(phantomSDName);
  sdMngr->AddNewDetector(phantomSD);
  layerLogical[2]->SetSensitiveDetector(phantomSD);

  //
  // primitive scorers
  //

  //                                        
  // Visualization attributes
  //
  //  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  phantomLogical->SetVisAttributes(simpleBoxVisAtt);

  G4VisAttributes * invisatt = new G4VisAttributes(G4Colour(.0,.0,.0));
  invisatt->SetVisibility(false);
  layerLogical[0]->SetVisAttributes(invisatt);
  layerLogical[1]->SetVisAttributes(invisatt);
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.8,.8,.8));
  visatt->SetVisibility(true);
  layerLogical[2]->SetVisAttributes(visatt);

}

void Test2DetectorConstruction::SetupDetectors()
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

